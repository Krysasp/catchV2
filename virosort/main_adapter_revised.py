# Obfuscated CATCHv2 script
# Original: main_adapter_revised.py
# Functions: ProbeMatch, reverse_complement, trim_adapters, bases_to_least_degenerate_iupac, consolidate_probe_sequence, match_probe_to_all_consensus, _match_single_probe, run_virosort_adapter_trim, main...

#!/usr/bin/env python3
"""
ViroSort Pipeline with Adapter Trimming and Multi-Cluster Probe Mapping

Processes probes by:
1. Removing adapter sequences from 5' and/or 3' ends
2. If adapter detection fails, trim fixed length from both ends
3. Map trimmed oligos to ALL consensus sequences with IUPAC/N tolerance
4. Consolidate matches to multiple clusters with "least degenerate" IUPAC sequence
5. Report start/end coordinates, cluster assignments, and taxonomy
"""

import argparse
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, NamedTuple
from multiprocessing import Pool, cpu_count

from scripts.parsers import (
    parse_fasta,
    parse_cluster_file,
    parse_probemap,
    parse_analysis_file
)
from scripts.taxonomy import load_csv_metadata


class ProbeMatch(NamedTuple):
    """Result of probe matching to multiple consensus sequences."""
    probe_id: str
    original_sequence: str
    consolidated_sequence: str
    adapters_removed: List[str]
    cluster_nums: List[int]
    cluster_coords: List[Tuple[int, int, int]]  # List of (cluster_num, start, end) tuples
    unmatched_bases: int
    taxonomy: Dict[str, str]
    organism_name: str
    num_sequences: int
    molecule_type: str
    reverse_complement: str
    mismatch_threshold: int


# IUPAC ambiguity codes and their base expansions
IUPAC_EXPANSION = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
    'n': ['A', 'C', 'G', 'T']
}

# Reverse mapping: specific bases to their IUPAC codes
BASE_TO_IUPAC = {
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

# Sort IUPAC codes by degeneracy (least degenerate first)
IUPAC_DEGENERACY = {
    'A': 1, 'C': 1, 'G': 1, 'T': 1,
    'S': 2, 'W': 2, 'K': 2, 'M': 2, 'R': 2, 'Y': 2,
    'B': 3, 'D': 3, 'H': 3, 'V': 3,
    'N': 4
}


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
        'N': 'N', 'n': 'N'
    }
    return ''.join(complement.get(base, base) for base in reversed(seq))


def trim_adapters(probe_seq: str, adapters: List[str], trim_len: int = 20) -> Tuple[str, List[str]]:
    """
    Remove adapter sequences from 5' and/or 3' ends of probe.
    
    Tries to match known adapters at both ends (including reverse complements).
    If no adapter found, falls back to trimming fixed length from both ends.
    
    Args:
        probe_seq: probe sequence
        adapters: list of adapter sequences
        trim_len: number of bases to trim from each end if adapter not detected
    
    Returns:
        tuple of (trimmed_sequence, list_of_removed_adapters)
    """
    original_seq = probe_seq
    trimmed = probe_seq
    removed_adapters = []
    
    # Create adapter patterns (forward and reverse complement)
    adapter_patterns = []
    for adapter in adapters:
        adapter_patterns.append(adapter)
        adapter_patterns.append(reverse_complement(adapter))
    
    # Try trimming from 5' end
    adapter_found_5 = False
    for adapter in adapter_patterns:
        if trimmed.startswith(adapter):
            trimmed = trimmed[len(adapter):]
            removed_adapters.append(f"5'{adapter}'")
            adapter_found_5 = True
            break
    
    # Try trimming from 3' end (use current trimmed sequence)
    adapter_found_3 = False
    for adapter in adapter_patterns:
        if trimmed.endswith(adapter):
            trimmed = trimmed[:-len(adapter)]
            removed_adapters.append(f"3'{adapter}'")
            adapter_found_3 = True
            break
    
    # Fallback: if no adapters found, trim fixed length from both ends
    if not adapter_found_5 and not adapter_found_3:
        trimmed = original_seq
        
        # Trim 5' end by fixed length
        if len(trimmed) > trim_len * 2:
            trimmed = trimmed[trim_len:]
            removed_adapters.append(f"5'fixed-{trim_len}bp")
            
            # Trim 3' end by fixed length
            if len(trimmed) > trim_len:
                removed_adapters.append(f"3'fixed-{trim_len}bp")
                trimmed = trimmed[:-trim_len]
    
    return trimmed, removed_adapters


def bases_to_least_degenerate_iupac(bases: set) -> str:
    """
    Convert a set of possible bases to the least degenerate IUPAC code.
    
    Args:
        bases: set of possible bases at a position
    
    Returns:
        IUPAC code representing the set (preferring less degenerate codes)
    """
    if not bases:
        return 'N'
    
    # Convert to frozenset for dictionary lookup
    base_frozenset = frozenset(bases)
    
    # Try direct match first
    if base_frozenset in BASE_TO_IUPAC:
        return BASE_TO_IUPAC[base_frozenset]
    
    # Fallback to N
    return 'N'


def consolidate_probe_sequence(trimmed_seq: str, consensus_matches: List[Tuple[int, int, int, int]], 
                               consensus_dict: Dict[int, str]) -> str:
    """
    Consolidate probe sequence from multiple consensus matches.
    
    For each position, collect all bases from matching consensus sequences,
    then use the least degenerate IUPAC code that represents that set.
    
    Args:
        trimmed_seq: original trimmed probe sequence
        consensus_matches: list of (cluster_num, start, end, unmatched) tuples
        consensus_dict: dictionary of {cluster_num: consensus_sequence}
    
    Returns:
        consolidated probe sequence with appropriate IUPAC codes
    """
    if not consensus_matches:
        return trimmed_seq
    
    probe_len = len(trimmed_seq)
    
    # Collect all bases at each position from all matching consensus sequences
    position_bases = [set() for _ in range(probe_len)]
    
    for cluster_num, start, end, unmatched in consensus_matches:
        consensus_seq = consensus_dict.get(cluster_num, '')
        if not consensus_seq:
            continue
        
        # Determine the position offset (consensus_seq may be longer than probe)
        # The trimmed_seq was found at position 'start' in consensus_seq
        for i in range(probe_len):
            if start + i < len(consensus_seq):
                c_base = consensus_seq[start + i].upper()
                # Get all possible bases for this IUPAC code
                if c_base in IUPAC_EXPANSION:
                    position_bases[i].update(IUPAC_EXPANSION[c_base])
                else:
                    position_bases[i].add(c_base)
    
    # Convert each position to the least degenerate IUPAC code
    consolidated = []
    for bases in position_bases:
        if len(bases) == 1:
            # Single base - use it directly
            consolidated.append(list(bases)[0])
        else:
            # Multiple bases - find the least degenerate IUPAC code
            iupac_code = bases_to_least_degenerate_iupac(bases)
            consolidated.append(iupac_code)
    
    return ''.join(consolidated)


def match_probe_to_all_consensus(probe_seq: str, consensus_dict: Dict[int, str], 
                                 max_unmatched: int = 6,
                                 window_size: int = 10,
                                 use_sliding_window: bool = False) -> List[Tuple[int, int, int, int]]:
    """
    Match trimmed oligo to ALL consensus sequences with IUPAC/N tolerance.
    Uses seed-based approach by default, with optional full sliding window.
    
    Returns a list of matches across all consensus sequences.
    
    Args:
        probe_seq: trimmed probe sequence
        consensus_dict: dictionary of {cluster_num: consensus_sequence}
        max_unmatched: maximum tolerated unmatched bases
        window_size: number of consensus sequences to search
        use_sliding_window: if True, use full sliding window search (slower but more thorough)
    
    Returns:
        list of (cluster_num, start, end, unmatched_bases) tuples for all matches
    """
    probe_len = len(probe_seq)
    probe_upper = probe_seq.upper()
    
    # Extract seed from probe (20bp contiguous non-IUPAC stretch)
    seed_bases = []
    seed_indices = []
    for i, base in enumerate(probe_upper):
        if base in ['A', 'C', 'G', 'T']:
            seed_bases.append(base)
            seed_indices.append(i)
            if len(seed_bases) >= 20:
                break
    
    sorted_clusters = sorted(consensus_dict.keys())
    all_matches = []
    
    # Search all clusters if window_size is None, otherwise limit to window_size
    search_clusters = sorted_clusters if window_size is None else sorted_clusters[:window_size]
    
    for cluster_num in search_clusters:
        consensus_seq = consensus_dict[cluster_num].upper()
        consensus_len = len(consensus_seq)
        
        if probe_len > consensus_len:
            continue
        
        if use_sliding_window:
            # Full sliding window search - check every position in consensus
            for start_pos in range(consensus_len - probe_len + 1):
                unmatched = 0
                for i in range(probe_len):
                    p_base = probe_upper[i]
                    c_base = consensus_seq[start_pos + i]
                    
                    if p_base in IUPAC_EXPANSION:
                        if c_base not in IUPAC_EXPANSION[p_base]:
                            unmatched += 1
                    elif p_base != c_base:
                        unmatched += 1
                    
                    if unmatched > max_unmatched:
                        break
                
                if unmatched <= max_unmatched:
                    all_matches.append((cluster_num, start_pos, start_pos + probe_len, unmatched))
        elif len(seed_bases) >= 15:
            # Seed-based approach - extract 15-20bp seed from probe
            seed_seq = ''.join(seed_bases)
            seed_start_idx = seed_indices[0]
            
            # Find all occurrences of seed in consensus
            search_start = 0
            while True:
                pos = consensus_seq.find(seed_seq, search_start)
                if pos == -1:
                    break
                
                actual_start = pos - seed_start_idx
                if actual_start < 0:
                    search_start = pos + 1
                    continue
                
                actual_end = actual_start + probe_len
                if actual_end > consensus_len:
                    search_start = pos + 1
                    continue
                
                # Verify match with full IUPAC-aware comparison
                unmatched = 0
                for i in range(probe_len):
                    p_base = probe_upper[i]
                    c_base = consensus_seq[actual_start + i]
                    
                    if p_base in IUPAC_EXPANSION:
                        if c_base not in IUPAC_EXPANSION[p_base]:
                            unmatched += 1
                    elif p_base != c_base:
                        unmatched += 1
                    
                    if unmatched > max_unmatched:
                        break
                
                if unmatched <= max_unmatched:
                    all_matches.append((cluster_num, actual_start, actual_end, unmatched))
                
                search_start = pos + 1
        else:
            # Fallback to sliding window for very short or highly degenerate probes
            for start_pos in range(consensus_len - probe_len + 1):
                unmatched = 0
                for i in range(probe_len):
                    p_base = probe_upper[i]
                    c_base = consensus_seq[start_pos + i]
                    
                    if p_base in IUPAC_EXPANSION:
                        if c_base not in IUPAC_EXPANSION[p_base]:
                            unmatched += 1
                    elif p_base != c_base:
                        unmatched += 1
                    
                    if unmatched > max_unmatched:
                        break
                
                if unmatched <= max_unmatched:
                    all_matches.append((cluster_num, start_pos, start_pos + probe_len, unmatched))
    
    return all_matches


def _match_single_probe(args: Tuple) -> Tuple[str, str, List[Tuple[int, int, int, int]], int]:
    """
    Helper function for parallel processing. Matches a single probe with all mismatch thresholds.
    
    Args:
        args: Tuple of (probe_id, probe_seq, consensus_dict, adapters, adapter_trim_enabled, trim_len, skip_alignment, cluster_seqs)
    
    Returns:
        Tuple of (probe_id, trimmed_seq, matches, mismatch_threshold) or (probe_id, trimmed_seq, [], threshold) if no match
    """
    probe_id, probe_seq, consensus_dict, adapters, adapter_trim_enabled, trim_len, skip_alignment, cluster_seqs = args
    
    # Trim adapters if enabled
    if adapter_trim_enabled:
        trimmed_seq, _ = trim_adapters(probe_seq, adapters, trim_len)
    else:
        trimmed_seq = probe_seq
    
    # Try to extract cluster number from probe header
    cluster_num_from_header = None
    cluster_match = re.search(r'CONS\|(\d+)\|', probe_id)
    if cluster_match:
        cluster_num_from_header = int(cluster_match.group(1))
    
    if skip_alignment and cluster_num_from_header is not None:
        # Skip alignment - use cluster number from header directly
        if cluster_num_from_header in consensus_dict:
            # Create a synthetic match using the cluster number from header
            probe_len = len(trimmed_seq)
            # Use start position 0 as placeholder since we're not doing alignment
            matches = [(cluster_num_from_header, 0, probe_len, 0)]
            return (probe_id, trimmed_seq, matches, 0)
    
    # Try mismatch thresholds 6, 7, 8, 9 (alignment mode)
    for mismatch_threshold in [6, 7, 8, 9]:
        matches = match_probe_to_all_consensus(
            trimmed_seq,
            consensus_dict,
            mismatch_threshold,
            window_size=None,
            use_sliding_window=True  # Always use sliding window for thoroughness
        )
        if matches:
            return (probe_id, trimmed_seq, matches, mismatch_threshold)
    
    return (probe_id, trimmed_seq, [], 9)  # No match found
    """
    Match trimmed oligo to ALL consensus sequences with IUPAC/N tolerance.
    Uses seed-based approach by default, with optional full sliding window.
    
    Returns a list of matches across all consensus sequences.
    
    Args:
        probe_seq: trimmed probe sequence
        consensus_dict: dictionary of {cluster_num: consensus_sequence}
        max_unmatched: maximum tolerated unmatched bases
        window_size: number of consensus sequences to search
        use_sliding_window: if True, use full sliding window search (slower but more thorough)
    
    Returns:
        list of (cluster_num, start, end, unmatched_bases) tuples for all matches
    """
    probe_len = len(probe_seq)
    probe_upper = probe_seq.upper()
    
    # Extract seed from probe (20bp contiguous non-IUPAC stretch)
    seed_bases = []
    seed_indices = []
    for i, base in enumerate(probe_upper):
        if base in ['A', 'C', 'G', 'T']:
            seed_bases.append(base)
            seed_indices.append(i)
            if len(seed_bases) >= 20:
                break
    
    sorted_clusters = sorted(consensus_dict.keys())
    all_matches = []
    
    # Search all clusters if window_size is None, otherwise limit to window_size
    search_clusters = sorted_clusters if window_size is None else sorted_clusters[:window_size]
    
    for cluster_num in search_clusters:
        consensus_seq = consensus_dict[cluster_num].upper()
        consensus_len = len(consensus_seq)
        
        if probe_len > consensus_len:
            continue
        
        if use_sliding_window:
            # Full sliding window search - check every position in consensus
            for start_pos in range(consensus_len - probe_len + 1):
                unmatched = 0
                for i in range(probe_len):
                    p_base = probe_upper[i]
                    c_base = consensus_seq[start_pos + i]
                    
                    if p_base in IUPAC_EXPANSION:
                        if c_base not in IUPAC_EXPANSION[p_base]:
                            unmatched += 1
                    elif p_base != c_base:
                        unmatched += 1
                    
                    if unmatched > max_unmatched:
                        break
                
                if unmatched <= max_unmatched:
                    all_matches.append((cluster_num, start_pos, start_pos + probe_len, unmatched))
        elif len(seed_bases) >= 15:
            # Seed-based approach - extract 15-20bp seed from probe
            seed_seq = ''.join(seed_bases)
            seed_start_idx = seed_indices[0]
            
            # Find all occurrences of seed in consensus
            search_start = 0
            while True:
                pos = consensus_seq.find(seed_seq, search_start)
                if pos == -1:
                    break
                
                actual_start = pos - seed_start_idx
                if actual_start < 0:
                    search_start = pos + 1
                    continue
                
                actual_end = actual_start + probe_len
                if actual_end > consensus_len:
                    search_start = pos + 1
                    continue
                
                # Verify match with full IUPAC-aware comparison
                unmatched = 0
                for i in range(probe_len):
                    p_base = probe_upper[i]
                    c_base = consensus_seq[actual_start + i]
                    
                    if p_base in IUPAC_EXPANSION:
                        if c_base not in IUPAC_EXPANSION[p_base]:
                            unmatched += 1
                    elif p_base != c_base:
                        unmatched += 1
                    
                    if unmatched > max_unmatched:
                        break
                
                if unmatched <= max_unmatched:
                    all_matches.append((cluster_num, actual_start, actual_end, unmatched))
                
                search_start = pos + 1
        else:
            # Fallback to sliding window for very short or highly degenerate probes
            for start_pos in range(consensus_len - probe_len + 1):
                unmatched = 0
                for i in range(probe_len):
                    p_base = probe_upper[i]
                    c_base = consensus_seq[start_pos + i]
                    
                    if p_base in IUPAC_EXPANSION:
                        if c_base not in IUPAC_EXPANSION[p_base]:
                            unmatched += 1
                    elif p_base != c_base:
                        unmatched += 1
                    
                    if unmatched > max_unmatched:
                        break
                
                if unmatched <= max_unmatched:
                    all_matches.append((cluster_num, start_pos, start_pos + probe_len, unmatched))
    
    return all_matches


def run_virosort_adapter_trim(
    degapped_fasta_path: str,
    cluster_file_path: str,
    probemap_path: str,
    analysis_path: str,
    metadata_path: str,
    output_dir: str,
    adapters: List[str] = None,
    max_unmatched: int = 6,
    trim_len: int = 20,
    max_clusters: int = None,
    use_sliding_window: bool = False,
    skip_alignment: bool = False,
    oligo_fasta_path: str = None
) -> Dict:
    """
    Main ViroSort pipeline with adapter trimming and multi-cluster probe mapping.
    """
    if adapters is None:
        adapters = [
            'GACCATCTAGCGACCTCCAC',
            'AGGCCCTGGCTGCTGATATG',
            'GACCATCTAGCGACCTCCTC',
            'GACCTTTTGGGACAGCGGTG'
        ]
    
    print("=" * 60)
    print("ViroSort Pipeline with Adapter Trimming")
    print("=" * 60)
    print(f"\nInput files:")
    print(f"  Degapped consensus: {degapped_fasta_path}")
    print(f"  Cluster file: {cluster_file_path}")
    print(f"  Probe map: {probemap_path}")
    print(f"  Analysis file: {analysis_path}")
    print(f"  Metadata: {metadata_path}")
    print(f"\nParameters:")
    print(f"  Adapters: {adapters}")
    print(f"  Trim length (fallback): {trim_len}bp")
    print(f"  Max unmatched bases: {max_unmatched}")
    print()
    
    print("Parsing input files...")
    degapped = parse_fasta(degapped_fasta_path)
    cluster_seqs = parse_cluster_file(cluster_file_path)
    
    # Check if probemap exists
    probes = None
    if probemap_path and os.path.exists(probemap_path):
        probes = parse_probemap(probemap_path)
        print(f"  Loaded {len(probes)} probes")
    else:
        print(f"  No probe map provided (or file not found: {probemap_path})")
    
    # Parse 26-mer oligo FASTA file if provided
    oligo_info = {}  # Maps cluster_num -> list of (oligo_id, oligo_seq, additional_info)
    if oligo_fasta_path and os.path.exists(oligo_fasta_path):
        oligo_fasta = parse_fasta(oligo_fasta_path)
        for header, seq in oligo_fasta.items():
            # Parse header: e.g., "da4a936606_CONS|24|220|65.078" or "5bcafea0c7_CONS|6|95pct85pct_degap"
            parts = header.split('|')
            if len(parts) >= 3:
                try:
                    cluster_num = int(parts[1])
                    # Extract additional info
                    if len(parts) >= 4:
                        # 26-mer format: id_CONS|cluster|pos|sim
                        additional_info = f"{parts[2]}|{parts[3]}"
                    else:
                        # Original format: id_CONS|cluster|95pct85pct_degap
                        additional_info = parts[2]
                    # Extract original identifier (before _CONS|)
                    header_parts = parts[0].split('_CONS')
                    original_id = header_parts[0] if header_parts else parts[0]
                    if cluster_num not in oligo_info:
                        oligo_info[cluster_num] = []
                    oligo_info[cluster_num].append({
                        'oligo_id': original_id,
                        'sequence': seq,
                        'additional_info': additional_info
                    })
                except ValueError:
                    pass
        print(f"  Loaded {len(oligo_info)} cluster entries from oligo FASTA")
    else:
        print(f"  No oligo FASTA provided (or file not found: {oligo_fasta_path})")
    
    analysis = parse_analysis_file(analysis_path)
    metadata = load_csv_metadata(metadata_path)
    print(f"  Loaded {len(degapped)} degapped consensus sequences")
    print(f"  Loaded {len(cluster_seqs)} clusters")
    print()
    
    print("Extracting consensus sequences...")
    degapped_consensus = {}
    for header, seq in degapped.items():
        parts = header.split('|')
        if len(parts) >= 2:
            try:
                cluster_num = int(parts[1])
                if max_clusters is None or cluster_num <= max_clusters:
                    degapped_consensus[cluster_num] = seq
            except ValueError:
                pass
    print(f"  Extracted {len(degapped_consensus)} consensus sequences")
    print()
    
    print("Building cluster-accession mapping...")
    cluster_accessions = {}
    cluster_molecule_types = {}
    for cluster_num, seq_list in cluster_seqs.items():
        cluster_int = int(cluster_num) if isinstance(cluster_num, str) else cluster_num
        if max_clusters is None or cluster_int <= max_clusters:
            cluster_accessions[cluster_int] = [seq['accession'] for seq in seq_list]
            # Get molecule type from first accession found in metadata
            for acc in seq_list:
                acc_id = acc['accession']
                # Try exact match first
                if acc_id in metadata:
                    mol_type = metadata[acc_id].get('molecule_type', 'Unknown')
                    cluster_molecule_types[cluster_int] = mol_type
                    break
                # Try matching by extracting base accession (before | if present)
                base_acc = acc_id.split('|')[0] if '|' in acc_id else acc_id
                if base_acc in metadata:
                    mol_type = metadata[base_acc].get('molecule_type', 'Unknown')
                    cluster_molecule_types[cluster_int] = mol_type
                    break
            else:
                cluster_molecule_types[cluster_int] = 'Unknown'
    print(f"  Mapped {len(cluster_accessions)} clusters to accessions")
    print()
    
    probe_matches: List[ProbeMatch] = []
    unmatched_probes = []
    adapter_stats = defaultdict(int)
    
    if probes:
        # Check if adapter trimming is disabled (adapters=['None'])
        adapter_trim_enabled = adapters != ['None']
        
        print("Processing probes (multi-cluster mapping)...")
        print(f"  Adapter trimming: {'Enabled' if adapter_trim_enabled else 'Disabled (using full probe sequence)'}")
        print(f"  Matching mode: {'Skip alignment (cluster from header)' if skip_alignment else 'Full sliding window (parallel)'}")
        if skip_alignment:
            print(f"  Will extract cluster number from probe headers (CONS|XXX|pattern)")
        else:
            print(f"  Testing mismatch thresholds: 6, 7, 8, 9")
        print()
        
        # Prepare probe data for parallel processing
        # Also build mapping from original identifier prefix to 26-mer oligo info
        probe_data = []
        probe_to_oligo_map = {}  # Maps probe_id (full) -> oligo_info dict
        for probe_id, probe_info in probes.items():
            original_seq = probe_info['sequence']
            num_seqs = probe_info.get('num_sequences', 0)
            probe_data.append((probe_id, original_seq, num_seqs))
            
            # Extract cluster number from probe header and map to 26-mer oligo info
            cluster_match = re.search(r'CONS\|(\d+)\|', probe_id)
            if cluster_match and oligo_info:
                cluster_num = int(cluster_match.group(1))
                # Find matching oligo in oligo_info for this cluster
                # The probe_id prefix (before _CONS|) should match oligo_id
                probe_prefix = probe_id.split('_CONS|')[0]
                for oligo_entry in oligo_info.get(cluster_num, []):
                    if oligo_entry['oligo_id'] == probe_prefix:
                        probe_to_oligo_map[probe_id] = oligo_entry
                        break
        
        # Use parallel processing with all available CPU cores
        num_cpus = cpu_count()
        print(f"  Using {num_cpus} CPU cores for parallel processing...")
        
        # Prepare arguments for parallel processing
        parallel_args = [
            (probe_id, seq, degapped_consensus, adapters, adapter_trim_enabled, trim_len, skip_alignment, cluster_seqs)
            for probe_id, seq, num_seqs in probe_data
        ]
        
        # Process probes in parallel
        with Pool(num_cpus) as pool:
            results = pool.map(_match_single_probe, parallel_args)
        
        # Process results
        for probe_id, trimmed_seq, all_matches, final_mismatch in results:
            # Find original probe data
            original_seq = None
            num_seqs = 0
            for pid, seq, ns in probe_data:
                if pid == probe_id:
                    original_seq = seq
                    num_seqs = ns
                    break
            
            if all_matches:
                # Find the best match (lowest unmatched bases)
                best_match = min(all_matches, key=lambda x: (x[3], x[0]))
                best_cluster_num, best_start, best_end, best_unmatched = best_match
                
                # Consolidate sequence across all matching consensus sequences
                consolidated_seq = consolidate_probe_sequence(trimmed_seq, all_matches, degapped_consensus)
                
                # Get taxonomy from first matching cluster
                taxonomy = {}
                organism_name = 'Unknown'
                for acc in cluster_accessions.get(best_cluster_num, []):
                    if acc in metadata:
                        taxonomy = metadata[acc]
                        organism_name = taxonomy.get('species_display', '').strip() or 'Unknown'
                        break
                
                if not taxonomy:
                    taxonomy = {
                        'family': 'Unknown',
                        'genus': 'Unknown',
                        'species': 'Unknown',
                        'molecule_type': 'Unknown'
                    }
                    organism_name = 'Unknown'
                
                # Collect all unique cluster numbers
                unique_clusters = sorted(set([m[0] for m in all_matches]))
                
                # Format adapters_removed for output
                adapters_str = 'N/A (no trimming)' if not adapter_trim_enabled else 'None'
                
                # Build cluster coordinates list (cluster_num, start, end) for all matches
                cluster_coords = [(m[0], m[1], m[2]) for m in all_matches]
                
                # Get molecule type
                molecule_type = cluster_molecule_types.get(best_cluster_num, 'Unknown')
                
                # Generate reverse complement if molecule type is ssRNA(+) or ssRNA(-) or ssDNA variants
                rc_sequence = ''
                if molecule_type in ['ssRNA(+)', 'ssRNA(-)', 'ssDNA(+)', 'ssDNA(-)', 'ssDNA', 'RNA', 'ssRNA']:
                    rc_sequence = reverse_complement(consolidated_seq)
                
                probe_matches.append(ProbeMatch(
                    probe_id=probe_id,
                    original_sequence=original_seq,
                    consolidated_sequence=consolidated_seq,
                    adapters_removed=adapters_str,
                    cluster_nums=unique_clusters,
                    cluster_coords=cluster_coords,
                    unmatched_bases=best_unmatched,
                    taxonomy=taxonomy,
                    organism_name=organism_name,
                    num_sequences=num_seqs,
                    molecule_type=molecule_type,
                    reverse_complement=rc_sequence,
                    mismatch_threshold=final_mismatch
                ))
            else:
                adapters_str = 'N/A (no trimming)' if not adapter_trim_enabled else 'None'
                
                # Find trimmed sequence
                if adapter_trim_enabled:
                    trimmed_seq, _ = trim_adapters(original_seq, adapters, trim_len)
                else:
                    trimmed_seq = original_seq
                
                unmatched_probes.append({
                    'probe_id': probe_id,
                    'original_sequence': original_seq,
                    'trimmed_sequence': trimmed_seq,
                    'adapters_removed': adapters_str,
                    'final_mismatch': final_mismatch,
                    'num_sequences': num_seqs
                })
        
        print(f"  Successfully matched probes: {len(probe_matches)}")
        print(f"  Unmatched probes: {len(unmatched_probes)}")
        print()
    
    print("Creating output directory structure...")
    output_base = output_dir
    fasta_dir = os.path.join(output_base, 'fasta')
    summary_dir = os.path.join(output_base, 'summary')
    html_dir = os.path.join(output_base, 'html_report')
    
    os.makedirs(fasta_dir, exist_ok=True)
    os.makedirs(summary_dir, exist_ok=True)
    os.makedirs(html_dir, exist_ok=True)
    print(f"  Output base: {output_base}")
    print()
    
    print("Writing renamed consensus FASTA...")
    consensus_fasta_path = os.path.join(fasta_dir, 'renamed_consensus.fasta')
    with open(consensus_fasta_path, 'w') as f:
        for cluster_num in sorted(degapped_consensus.keys()):
            consensus_seq = degapped_consensus[cluster_num]
            
            taxonomy = {}
            for acc in cluster_accessions.get(cluster_num, []):
                if acc in metadata:
                    taxonomy = metadata[acc]
                    break
            
            if not taxonomy:
                taxonomy = {'family': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown'}
            
            coverage = '0.000'
            for genome_id, stats in analysis.items():
                if f'genome {cluster_num}' in genome_id:
                    coverage = f"{stats.get('avg_coverage_over_unambig', 0):.3f}"
                    break
            
            family = taxonomy.get('family', 'Unknown').replace(' ', '_')
            genus = taxonomy.get('genus', 'Unknown').replace(' ', '_')
            species = taxonomy.get('species', 'Unknown').replace(' ', '_')
            molecule_type = cluster_molecule_types.get(cluster_num, 'Unknown')
            
            header = f"C{int(cluster_num):04d}|{family}|{genus}|{species}|{coverage}|{molecule_type}"
            f.write(f'>{header}\n{consensus_seq}\n')
    print(f"  Written: {consensus_fasta_path}")
    
    # Only write probes if they were provided
    if probes and probe_matches:
        print("Writing renamed probes FASTA...")
        probes_fasta_path = os.path.join(fasta_dir, 'renamed_probes.fasta')
        with open(probes_fasta_path, 'w') as f:
            for match in probe_matches:
                taxonomy = match.taxonomy
                family = taxonomy.get('family', 'Unknown').replace(' ', '_')
                genus = taxonomy.get('genus', 'Unknown').replace(' ', '_')
                species = taxonomy.get('species', 'Unknown').replace(' ', '_')
                molecule_type = match.molecule_type
                coverage = '0.000'
                for genome_id, stats in analysis.items():
                    if f'genome {match.cluster_nums[0]}' in genome_id:
                        coverage = f"{stats.get('avg_coverage_over_unambig', 0):.3f}"
                        break
                
                # Extract probe identifier prefix and additional info from probe_id
                # Header format: "4d97fd76cd_CONS|24|95pct85pct_degap"
                # Split by _CONS to get the identifier, then parse the rest
                header_parts = match.probe_id.split('_CONS')
                probe_prefix = header_parts[0] if header_parts else match.probe_id
                # The rest after _CONS is "|24|95pct85pct_degap", split by | to get cluster and additional_info
                rest_parts = header_parts[1].split('|') if len(header_parts) > 1 else []
                additional_info = rest_parts[2] if len(rest_parts) > 2 else ''
                
                # Use the original probe sequence from probemap (not the 26-mer window)
                probe_seq = match.original_sequence
                
                # Format cluster numbers as underscore-separated list
                cluster_str = '_'.join([f'C{int(c):04d}' for c in match.cluster_nums])
                header = f"{probe_prefix}|{cluster_str}|{additional_info}|{family}|{genus}|{species}|{match.num_sequences}|{coverage}|{molecule_type}"
                f.write(f'>{header}\n{probe_seq}\n')
                # Write reverse complement for ssRNA/ssDNA probes
                if match.reverse_complement:
                    rc_header = f"{probe_prefix}_RC|{cluster_str}|{additional_info}|{family}|{genus}|{species}|{match.num_sequences}|{coverage}|{molecule_type}"
                    f.write(f'>{rc_header}\n{match.reverse_complement}\n')
        print(f"  Written: {probes_fasta_path} ({len(probe_matches)} probes)")
        
        # Write 26-mer oligo FASTA with original identifiers
        print("Writing 26-mer oligo FASTA with original identifiers...")
        oligo_fasta_path = os.path.join(fasta_dir, 'renamed_26mer_oligos.fasta')
        with open(oligo_fasta_path, 'w') as f:
            # Iterate through oligo_info directly to get original 26-mer entries
            for cluster_num in sorted(oligo_info.keys()):
                for oligo_entry in oligo_info[cluster_num]:
                    oligo_id = oligo_entry['oligo_id']
                    oligo_seq = oligo_entry['sequence']
                    additional_info = oligo_entry['additional_info']
                    
                    # Get taxonomy from metadata using cluster accessions
                    taxonomy = {}
                    molecule_type = cluster_molecule_types.get(cluster_num, 'Unknown')
                    for acc in cluster_accessions.get(cluster_num, []):
                        if acc in metadata:
                            taxonomy = metadata[acc]
                            break
                    
                    if not taxonomy:
                        taxonomy = {'family': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown'}
                    
                    family = taxonomy.get('family', 'Unknown').replace(' ', '_')
                    genus = taxonomy.get('genus', 'Unknown').replace(' ', '_')
                    species = taxonomy.get('species', 'Unknown').replace(' ', '_')
                    
                    coverage = '0.000'
                    for genome_id, stats in analysis.items():
                        if f'genome {cluster_num}' in genome_id:
                            coverage = f"{stats.get('avg_coverage_over_unambig', 0):.3f}"
                            break
                    
                    # Format cluster number
                    cluster_str = f'C{int(cluster_num):04d}'
                    num_seqs = len(cluster_accessions.get(cluster_num, []))
                    
                    header = f"{oligo_id}|{cluster_str}|{additional_info}|{family}|{genus}|{species}|{num_seqs}|{coverage}|{molecule_type}"
                    f.write(f'>{header}\n{oligo_seq}\n')
                    # Write reverse complement for ssRNA/ssDNA
                    if molecule_type in ['ssRNA(+)', 'ssRNA(-)', 'ssDNA(+)', 'ssDNA(-)', 'ssDNA', 'RNA', 'ssRNA']:
                        rc_sequence = reverse_complement(oligo_seq)
                        rc_header = f"{oligo_id}_RC|{cluster_str}|{additional_info}|{family}|{genus}|{species}|{num_seqs}|{coverage}|{molecule_type}"
                        f.write(f'>{rc_header}\n{rc_sequence}\n')
        print(f"  Written: {oligo_fasta_path}")
        
        print("Writing cluster summary TSV...")
        summary_tsv_path = os.path.join(summary_dir, 'cluster_summary.tsv')
        with open(summary_tsv_path, 'w') as f:
            f.write('Cluster\tFamily\tGenus\tSpecies\tMolecule_Type\tCoverage\tNum_Probes_Matched\tConsensus_Length\n')
            for cluster_num in sorted(degapped_consensus.keys()):
                taxonomy = {}
                for acc in cluster_accessions.get(cluster_num, []):
                    if acc in metadata:
                        taxonomy = metadata[acc]
                        break
                
                if not taxonomy:
                    taxonomy = {'family': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown', 'molecule_type': 'Unknown'}
                
                coverage = '0.000'
                for genome_id, stats in analysis.items():
                    if f'genome {cluster_num}' in genome_id:
                        coverage = f"{stats.get('avg_coverage_over_unambig', 0):.3f}"
                        break
                
                num_probes = sum(1 for m in probe_matches if cluster_num in m.cluster_nums)
                consensus_length = len(degapped_consensus[cluster_num])
                molecule_type = cluster_molecule_types.get(cluster_num, 'Unknown')
                
                f.write(f"{cluster_num}\t{taxonomy.get('family', 'Unknown')}\t"
                       f"{taxonomy.get('genus', 'Unknown')}\t{taxonomy.get('species', 'Unknown')}\t"
                       f"{molecule_type}\t{coverage}\t{num_probes}\t{consensus_length}\n")
        print(f"  Written: {summary_tsv_path}")
        
        print("Writing probe mapping Excel...")
        excel_path = os.path.join(summary_dir, 'probe_mapping.xlsx')
        
        try:
            from openpyxl import Workbook
            wb = Workbook()
            ws = wb.active
            ws.title = "Probe_Mapping"
            
            headers = [
                'Probe_ID', 'Original_Sequence', 'Consolidated_Sequence', 'Adapters_Removed',
                'Cluster_Nums', 'Start_Pos', 'End_Pos', 'Unmatched_Bases',
                'Family', 'Genus', 'Species', 'Molecule_Type', 'Coverage', 'Num_Sequences'
            ]
            ws.append(headers)
            
            for match in probe_matches:
                adapters_str = match.adapters_removed if match.adapters_removed else 'None'
                cluster_str = '_'.join([f'C{int(c):04d}' for c in match.cluster_nums])
                ws.append([
                    match.probe_id,
                    match.original_sequence,
                    match.consolidated_sequence,
                    adapters_str,
                    cluster_str,
                    match.start,
                    match.end,
                    match.unmatched_bases,
                    match.taxonomy.get('family', 'Unknown'),
                    match.taxonomy.get('genus', 'Unknown'),
                    match.taxonomy.get('species', 'Unknown'),
                    cluster_molecule_types.get(match.cluster_nums[0], 'Unknown'),
                    coverage,
                    match.num_sequences
                ])
            
            for col in ws.columns:
                max_len = max(len(str(cell.value)) for cell in col) if col else 10
                col[0].column_letter = chr(65 + col[0].column - 1)
                ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 30)
            
            wb.save(excel_path)
            print(f"  Written: {excel_path}")
        except ImportError:
            print(f"  openpyxl not available, writing TSV fallback...")
            excel_tsv_path = os.path.join(summary_dir, 'probe_mapping.tsv')
            with open(excel_tsv_path, 'w') as f:
                f.write('Probe_ID\tOriginal_Sequence\tConsolidated_Sequence\tReverse_Complement\t')
                f.write('Adapters_Removed\tCluster_Nums\tCluster_Coords\tMolecule_Type\tMismatch_Threshold\t')
                f.write('Unmatched_Bases\tFamily\tGenus\tSpecies\tOrganism_Name\tCoverage\tNum_Sequences\n')
                for match in probe_matches:
                    adapters_str = match.adapters_removed if match.adapters_removed else 'None'
                    cluster_str = ', '.join([f'C{int(c):04d}' for c in match.cluster_nums])
                    coord_str = ', '.join([f'{coord[1]}-{coord[2]}' for coord in match.cluster_coords])
                    f.write(f"{match.probe_id}\t{match.original_sequence}\t{match.consolidated_sequence}\t")
                    f.write(f"{match.reverse_complement}\t{adapters_str}\t{cluster_str}\t{coord_str}\t")
                    f.write(f"{match.molecule_type}\t{match.mismatch_threshold}\t{match.unmatched_bases}\t")
                    f.write(f"{match.taxonomy.get('family', 'Unknown')}\t{match.taxonomy.get('genus', 'Unknown')}\t")
                    f.write(f"{match.taxonomy.get('species', 'Unknown')}\t{match.organism_name}\t0.000\t{match.num_sequences}\n")
            print(f"  Written (TSV fallback): {excel_tsv_path}")
        
        # Write summary CSV with grouped statistics
        print("Writing summary CSV...")
        summary_csv_path = os.path.join(summary_dir, 'summary.csv')
        with open(summary_csv_path, 'w') as f:
            f.write('Group_Type,Group_Name,Molecule_Type,Total_Consensus,Unique_Clusters,'
                   f'Total_Probes,Probes_Matched,Avg_Mismatch_Threshold\n')
            
            # Group by molecule type
            mol_type_groups = defaultdict(lambda: {
                'consensus': set(), 'clusters': set(), 'probes': set(), 
                'mismatch_sum': 0, 'count': 0
            })
            
            for match in probe_matches:
                mol_type = match.molecule_type
                mol_type_groups[mol_type]['consensus'].update(match.cluster_nums)
                mol_type_groups[mol_type]['clusters'].update(match.cluster_nums)
                mol_type_groups[mol_type]['probes'].add(match.probe_id)
                mol_type_groups[mol_type]['mismatch_sum'] += match.mismatch_threshold
                mol_type_groups[mol_type]['count'] += 1
            
            for mol_type in sorted(mol_type_groups.keys()):
                stats = mol_type_groups[mol_type]
                avg_mismatch = stats['mismatch_sum'] / stats['count'] if stats['count'] > 0 else 0
                f.write(f"Cluster,{mol_type},{mol_type},{len(stats['consensus'])},"
                       f"{len(stats['clusters'])},{len(stats['probes'])},{stats['count']},"
                       f"{avg_mismatch:.2f}\n")
        
        print(f"  Written: {summary_csv_path}")
        
        print("Writing HTML visualization...")
        html_path = os.path.join(html_dir, 'probe_mapping_report.html')
        
        with open(html_path, 'w') as f:
            f.write('<!DOCTYPE html>\n')
            f.write('<html><head>\n')
            f.write('<title>ViroSort Probe Mapping Report</title>\n')
            f.write('<style>\n')
            f.write('body { font-family: Arial, sans-serif; margin: 20px; }\n')
            f.write('h1, h2 { color: #333; }\n')
            f.write('table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n')
            f.write('th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n')
            f.write('th { background-color: #4CAF50; color: white; }\n')
            f.write('tr:nth-child(even) { background-color: #f2f2f2; }\n')
            f.write('.sequence { font-family: monospace; word-break: break-all; }\n')
            f.write('.summary-box { background: #e7f3fe; padding: 15px; border-radius: 5px; margin: 20px 0; }\n')
            f.write('</style>\n')
            f.write('</head><body>\n')
            
            f.write('<h1>ViroSort Probe Mapping Report</h1>\n')
            
            f.write('<div class="summary-box">\n')
            f.write(f'<h2>Summary</h2>\n')
            f.write(f'<p><strong>Total Clusters:</strong> {len(degapped_consensus)}</p>\n')
            f.write(f'<p><strong>Total Probes Processed:</strong> {len(probes)}</p>\n')
            f.write(f'<p><strong>Successfully Matched Probes:</strong> {len(probe_matches)}</p>\n')
            f.write(f'<p><strong>Unmatched Probes:</strong> {len(unmatched_probes)}</p>\n')
            f.write(f'<p><strong>Mismatch Thresholds Tested:</strong> 6, 7, 8, 9</p>\n')
            f.write('</div>\n')
            
            f.write('<h2>Cluster Summary</h2>\n')
            f.write('<table>\n')
            f.write('<tr><th>Cluster</th><th>Family</th><th>Genus</th><th>Species</th><th>Probes Matched</th></tr>\n')
            for cluster_num in sorted(degapped_consensus.keys()):
                taxonomy = {}
                for acc in cluster_accessions.get(cluster_num, []):
                    if acc in metadata:
                        taxonomy = metadata[acc]
                        break
                num_probes = sum(1 for m in probe_matches if cluster_num in m.cluster_nums)
                f.write(f'<tr><td>C{int(cluster_num):04d}</td><td>{taxonomy.get("family", "Unknown")}</td>'
                       f'<td>{taxonomy.get("genus", "Unknown")}</td><td>{taxonomy.get("species", "Unknown")}</td>'
                       f'<td>{num_probes}</td></tr>\n')
            f.write('</table>\n')
            
            # Write summary statistics (grouped by cluster and molecule type)
            f.write('<h2>Summary Statistics</h2>\n')
            f.write('<table>\n')
            f.write('<tr><th>Molecule_Type</th><th>Total_Consensus</th><th>Total_Probes</th>'
                   f'<th>Unique_Clusters</th><th>Probes_Matched</th></tr>\n')
            
            mol_type_stats = defaultdict(lambda: {'consensus': set(), 'probes': set(), 'clusters': set(), 'matched': 0})
            for match in probe_matches:
                mol_type = match.molecule_type
                mol_type_stats[mol_type]['consensus'].update(match.cluster_nums)
                mol_type_stats[mol_type]['probes'].add(match.probe_id)
                mol_type_stats[mol_type]['clusters'].update(match.cluster_nums)
                mol_type_stats[mol_type]['matched'] += 1
            
            for mol_type in sorted(mol_type_stats.keys()):
                stats = mol_type_stats[mol_type]
                f.write(f'<tr><td>{mol_type}</td><td>{len(stats["consensus"])}</td>'
                       f'<td>{len(stats["probes"])}</td><td>{len(stats["clusters"])}</td>'
                       f'<td>{stats["matched"]}</td></tr>\n')
            f.write('</table>\n')
            
            # Write mismatch threshold distribution
            f.write('<h2>Mismatch Threshold Distribution</h2>\n')
            f.write('<table>\n')
            f.write('<tr><th>Mismatches</th><th>Probes_Matched</th><th>Percentage</th></tr>\n')
            mismatch_dist = defaultdict(int)
            for match in probe_matches:
                mismatch_dist[match.mismatch_threshold] += 1
            total_matched = len(probe_matches)
            for mismatches in sorted(mismatch_dist.keys()):
                percentage = (mismatch_dist[mismatches] / total_matched * 100) if total_matched > 0 else 0
                f.write(f'<tr><td>{mismatches}</td><td>{mismatch_dist[mismatches]}</td>'
                       f'<td>{percentage:.1f}%</td></tr>\n')
            f.write('</table>\n')
            
            f.write('<h2>Probe Mapping Details (First 100)</h2>\n')
            f.write('<table>\n')
            f.write('<tr><th>Probe ID</th><th>Clusters</th><th>Cluster Coordinates</th>'
                   f'<th>Molecule_Type</th><th>Mismatches</th><th>Consolidated Seq</th></tr>\n')
            for match in sorted(probe_matches, key=lambda x: (x.cluster_nums[0], x.cluster_coords[0][1] if x.cluster_coords else 0))[:100]:
                cluster_str = ', '.join([f'C{int(c):04d}' for c in match.cluster_nums])
                coord_str = ', '.join([f'{coord[1]}-{coord[2]}' for coord in match.cluster_coords])
                f.write(f'<tr><td>{match.probe_id}</td><td>{cluster_str}</td>'
                       f'<td>{coord_str}</td><td>{match.molecule_type}</td>'
                       f'<td>{match.mismatch_threshold}</td>'
                       f'<td class="sequence">{match.consolidated_sequence}</td></tr>\n')
            f.write('</table>\n')
            
            f.write('</body></html>\n')
        print(f"  Written: {html_path}")
    
    print()
    print("=" * 60)
    print("Pipeline complete!")
    print("=" * 60)
    
    return {
        'probe_matches': probe_matches,
        'unmatched_probes': unmatched_probes,
        'adapter_stats': adapter_stats,
        'output_dir': output_base
    }


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description='ViroSort pipeline with adapter trimming and multi-cluster probe mapping'
    )
    
    parser.add_argument('--degapped', required=True,
                       help='Path to degapped consensus FASTA file')
    parser.add_argument('--cluster', required=True,
                       help='Path to cluster file')
    parser.add_argument('--probes', default=None,
                       help='Path to probe map TSV file (optional)')
    parser.add_argument('--analysis', required=True,
                       help='Path to analysis file')
    parser.add_argument('--metadata', required=True,
                       help='Path to metadata CSV file')
    parser.add_argument('--output', required=True,
                       help='Output directory path')
    parser.add_argument('--adapters', nargs='+', default=None,
                       help='Adapter sequences to trim (default: 4 known adapters)')
    parser.add_argument('--trim-len', type=int, default=20,
                       help='Number of bases to trim from each end if adapter not detected (default: 20)')
    parser.add_argument('--max-unmatched', type=int, default=6,
                       help='Maximum unmatched bases allowed (default: 6)')
    parser.add_argument('--max-clusters', type=int, default=None,
                       help='Maximum number of clusters to process')
    parser.add_argument('--sliding-window', action='store_true', default=False,
                       help='Use full sliding window search (slower but more thorough)')
    parser.add_argument('--skip-alignment', action='store_true', default=False,
                       help='Skip sequence alignment; use cluster number from probe header (CONS|XXX|) directly')
    parser.add_argument('--oligos', default=None,
                       help='Path to 26-mer oligo FASTA file (optional)')
    
    args = parser.parse_args()
    
    result = run_virosort_adapter_trim(
        degapped_fasta_path=args.degapped,
        cluster_file_path=args.cluster,
        probemap_path=args.probes,
        analysis_path=args.analysis,
        metadata_path=args.metadata,
        output_dir=args.output,
        adapters=args.adapters,
        max_unmatched=args.max_unmatched,
        trim_len=args.trim_len,
        max_clusters=args.max_clusters,
        use_sliding_window=args.sliding_window,
        skip_alignment=args.skip_alignment,
        oligo_fasta_path=args.oligos
    )
    
    return result


if __name__ == '__main__':
    main()
