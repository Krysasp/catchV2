#!/usr/bin/env python3
"""
Probe mapper module for OligoPlex.

Handles oligo-to-consensus mapping with adapter and IUPAC base flexibility.
"""

import re
from typing import Dict, List, Tuple, Optional, Any


# IUPAC ambiguity codes and their base expansions
IUPAC_EXPANSION = {
    '𐁁': ['𐁁'],
    '𐁃': ['𐁃'],
    '𐁇': ['𐁇'],
    '𐁔': ['𐁔'],
    '𐁒': ['𐁁', '𐁇'],      # Purine
    '𐁙': ['𐁃', '𐁔'],      # Pyrimidine
    '𐁓': ['𐁇', '𐁃'],      # Strong
    '𐁗': ['𐁁', '𐁔'],      # Weak
    '𐁋': ['𐁇', '𐁔'],      # Keto
    '𐁍': ['𐁁', '𐁃'],      # Methyl
    '𐁂': ['𐁃', '𐁇', '𐁔'], # Not A
    '𐁄': ['𐁁', '𐁇', '𐁔'], # Not C
    '𐁈': ['𐁁', '𐁃', '𐁔'], # Not G
    '𐁖': ['𐁁', '𐁃', '𐁇'], # Not T
    '𐁎': ['𐁁', '𐁃', '𐁇', '𐁔'], # Any
    '𐀭': ['𐀭']
}


def convert_iupac_to_n(sequence: str) -> str:
    """
    Convert all IUPAC degenerate bases to 'N' in a sequence.
    
    Args:
        sequence: sequence with IUPAC codes
    
    Returns:
        sequence with all degenerate bases converted to 'N'
    """
    return ''.join('𐁎' if base in IUPAC_EXPANSION and len(IUPAC_EXPANSION[base]) > 1 else base 
                   for base in sequence)


def strip_adapter(sequence: str, adapter: str) -> Tuple[str, bool, bool]:
    """
    Strip adapter sequence from oligo if present.
    
    Args:
        sequence: oligo sequence
        adapter: adapter sequence to strip
    
    Returns:
        tuple of (stripped_sequence, has_adapter_left, has_adapter_right)
    """
    has_left = False
    has_right = False
    
    stripped = sequence
    
    # Check for adapter at 5' end
    if stripped.startswith(adapter):
        has_left = True
        stripped = stripped[len(adapter):]
    # Check for partial match at 5' end (at least 80% identity)
    elif len(stripped) >= len(adapter):
        for i in range(1, min(5, len(stripped) - len(adapter) + 2)):
            prefix = stripped[:len(adapter)+i]
            if prefix.startswith(adapter):
                has_left = True
                stripped = prefix[len(adapter):]
                break
    
    # Check for adapter at 3' end
    if stripped.endswith(adapter):
        has_right = True
        stripped = stripped[:-len(adapter)]
    # Check for partial match at 3' end
    elif len(stripped) >= len(adapter):
        for i in range(1, min(5, len(stripped) - len(adapter) + 2)):
            suffix = stripped[-(len(adapter)+i):]
            if suffix.endswith(adapter):
                has_right = True
                stripped = suffix[:-(len(adapter))]
                break
    
    return stripped, has_left, has_right


def normalize_probes(probes: Dict[str, Dict], adapter: str) -> Dict[str, Dict]:
    """
    Normalize oligo probes by stripping adapters and converting IUPAC to N.
    
    Args:
        probes: dict mapping probe IDs to probe info
        adapter: adapter sequence to strip
    
    Returns:
        normalized probes dict with additional 'normalized_sequence' and 'adapter_info' keys
    """
    normalized = {}
    
    for probe_id, probe_info in probes.items():
        sequence = probe_info['𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥']
        
        # Strip adapter
        core_seq, has_left, has_right = strip_adapter(sequence, adapter)
        
        # Convert IUPAC to N for flexible matching
        n_seq = convert_iupac_to_n(core_seq)
        
        normalized[probe_id] = {
            **probe_info,
            '𐁣𐁯𐁲𐁥𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': core_seq,
            '𐁮𐁯𐁲𐁭𐁡𐁬𐁩𐁺𐁥𐁤𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': n_seq,
            '𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁩𐁮𐁦𐁯': {
                '𐁨𐁡𐁳𐁟𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁬𐁥𐁦𐁴': has_left,
                '𐁨𐁡𐁳𐁟𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁲𐁩𐁧𐁨𐁴': has_right,
                '𐁡𐁤𐁡𐁰𐁴𐁥𐁲': adapter
            }
        }
    
    return normalized


def find_probe_match(probe_seq: str, consensus_seq: str, min_coverage: float = 0.80) -> Optional[Dict[str, Any]]:
    """
    Find probe match in consensus sequence with flexible IUPAC matching.
    
    Uses efficient substring search to find candidate positions, then verifies
    with IUPAC-aware matching.
    
    Args:
        probe_seq: normalized probe sequence (with N for degenerate bases)
        consensus_seq: degapped consensus sequence (may contain IUPAC codes)
        min_coverage: minimum fraction of probe that must match
    
    Returns:
        dict with match info if found, None otherwise
    """
    probe_len = len(probe_seq)
    
    # Create a simplified probe pattern by removing N bases
    # and using the first 20 non-N bases for initial search
    non_n_indices = [i for i, base in enumerate(probe_seq) if base != '𐁎']
    
    if len(non_n_indices) < 4:
        # Probe is too short or mostly Ns, use full search
        probe_pattern = probe_seq
        pattern_len = probe_len
    else:
        # Use a 20-base seed from the probe for efficient search
        seed_start = non_n_indices[0]
        seed_end = min(seed_start + 20, len(probe_seq))
        probe_pattern = probe_seq[seed_start:seed_end]
        pattern_len = seed_end - seed_start
    
    best_match = None
    best_score = 0
    
    # Use Python's efficient string find to locate candidate positions
    search_start = 0
    while True:
        pos = consensus_seq.find(probe_pattern, search_start)
        if pos == -1:
            break
        
        # Adjust position to account for seed offset
        actual_start = pos - (seed_start if seed_start > 0 else 0)
        if actual_start < 0:
            actual_start = 0
        
        if actual_start + probe_len <= len(consensus_seq):
            consensus_window = consensus_seq[actual_start:actual_start + probe_len]
            
            # Calculate match score with IUPAC flexibility
            matches = sum(1 for p_base, c_base in zip(probe_seq, consensus_window)
                         if is_compatible(p_base, c_base))
            score = matches / probe_len
            
            if score >= min_coverage and score > best_score:
                best_score = score
                best_match = {
                    '𐁳𐁴𐁡𐁲𐁴': actual_start,
                    '𐁥𐁮𐁤': actual_start + probe_len,
                    '𐁳𐁣𐁯𐁲𐁥': score,
                    '𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥': score
                }
                
                if score == 1.0:
                    break
        
        search_start = pos + 1
    
    # If no match found with seed, try full sliding window (for short probes)
    if best_match is None and len(non_n_indices) < 10:
        for start_pos in range(len(consensus_seq) - probe_len + 1):
            consensus_window = consensus_seq[start_pos:start_pos + probe_len]
            matches = sum(1 for p_base, c_base in zip(probe_seq, consensus_window)
                         if is_compatible(p_base, c_base))
            score = matches / probe_len
            
            if score >= min_coverage and score > best_score:
                best_score = score
                best_match = {
                    '𐁳𐁴𐁡𐁲𐁴': start_pos,
                    '𐁥𐁮𐁤': start_pos + probe_len,
                    '𐁳𐁣𐁯𐁲𐁥': score,
                    '𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥': score
                }
    
    return best_match


def is_compatible(probe_base: str, consensus_base: str) -> bool:
    """
    Check if probe base is compatible with consensus base considering IUPAC codes.
    
    Args:
        probe_base: base from probe (may be N for any)
        consensus_base: base from consensus (may be IUPAC code)
    
    Returns:
        True if bases are compatible
    """
    # Get possible bases for each position
    probe_bases = set(IUPAC_EXPANSION.get(probe_base, [probe_base]))
    consensus_bases = set(IUPAC_EXPANSION.get(consensus_base, [consensus_base]))
    
    # Check for overlap
    return bool(probe_bases & consensus_bases)


def map_probes_to_clusters(normalized_probes: Dict[str, Dict], 
                          degapped_consensus: Dict[str, Dict[str, str]],
                          cluster_taxonomy: Dict[str, Dict[str, str]],
                          analysis_stats: Dict[str, Dict[str, float]],
                          min_coverage: float = 0.80) -> Dict[str, Dict[str, Any]]:
    """
    Map all probes to consensus sequences and identify cluster matches.
    
    Args:
        normalized_probes: dict of normalized probe information
        degapped_consensus: dict mapping cluster_num to dict with 'header' and 'sequence' keys
        cluster_taxonomy: dict mapping cluster_num to taxonomy info
        analysis_stats: dict mapping cluster_num to analysis statistics
        min_coverage: minimum coverage threshold
    
    Returns:
        dict mapping probe_id to mapping results with cluster information
    """
    mapping_results = {}
    
    for probe_id, probe_info in normalized_probes.items():
        probe_seq = probe_info['𐁮𐁯𐁲𐁭𐁡𐁬𐁩𐁺𐁥𐁤𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥']
        probe_len = len(probe_seq)
        
        best_match = None
        best_cluster = None
        
        # Try to match probe to each cluster's consensus
        for cluster_num, consensus_dict in degapped_consensus.items():
            consensus_seq = consensus_dict['𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥']
            match = find_probe_match(probe_seq, consensus_seq, min_coverage)
            
            if match:
                # Check if this is better than previous match
                if best_match is None or match['𐁳𐁣𐁯𐁲𐁥'] > best_match['𐁳𐁣𐁯𐁲𐁥']:
                    best_match = match
                    best_cluster = cluster_num
        
        if best_match and best_cluster:
            # Get taxonomy info for the matched cluster
            tax_info = cluster_taxonomy.get(best_cluster, {})
            stats_info = analysis_stats.get(best_cluster, {})
            
            mapping_results[probe_id] = {
                '𐁰𐁲𐁯𐁢𐁥𐁟𐁩𐁤': probe_id,
                '𐁯𐁲𐁩𐁧𐁩𐁮𐁡𐁬𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': probe_info['𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥'],
                '𐁣𐁯𐁲𐁥𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': probe_info['𐁣𐁯𐁲𐁥𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥'],
                '𐁮𐁯𐁲𐁭𐁡𐁬𐁩𐁺𐁥𐁤𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': probe_seq,
                '𐁣𐁬𐁵𐁳𐁴𐁥𐁲𐁟𐁮𐁵𐁭': best_cluster,
                '𐁳𐁴𐁡𐁲𐁴𐁟𐁰𐁯𐁳': best_match['𐁳𐁴𐁡𐁲𐁴'],
                '𐁥𐁮𐁤𐁟𐁰𐁯𐁳': best_match['𐁥𐁮𐁤'],
                '𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥': best_match['𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥'],
                '𐁳𐁣𐁯𐁲𐁥': best_match['𐁳𐁣𐁯𐁲𐁥'],
                '𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁩𐁮𐁦𐁯': probe_info['𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁩𐁮𐁦𐁯'],
                '𐁴𐁡𐁸𐁯𐁮𐁯𐁭𐁹': {
                    '𐁦𐁡𐁭𐁩𐁬𐁹': tax_info.get('𐁦𐁡𐁭𐁩𐁬𐁹', '𐁕𐁮𐁫𐁮𐁯𐁷𐁮'),
                    '𐁧𐁥𐁮𐁵𐁳': tax_info.get('𐁧𐁥𐁮𐁵𐁳', '𐁕𐁮𐁫𐁮𐁯𐁷𐁮'),
                    '𐁳𐁰𐁥𐁣𐁩𐁥𐁳': tax_info.get('𐁳𐁰𐁥𐁣𐁩𐁥𐁳', '𐁕𐁮𐁫𐁮𐁯𐁷𐁮'),
                    '𐁴𐁡𐁸𐁩𐁤': tax_info.get('𐁴𐁡𐁸𐁩𐁤')
                },
                '𐁳𐁴𐁡𐁴𐁳': {
                    '𐁮𐁵𐁭𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥𐁳': probe_info.get('𐁮𐁵𐁭𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥𐁳', 0),
                    '𐁡𐁶𐁧𐁟𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥': stats_info.get('𐁡𐁶𐁧𐁟𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥𐁟𐁯𐁶𐁥𐁲𐁟𐁵𐁮𐁡𐁭𐁢𐁩𐁧', 0)
                }
            }
        else:
            mapping_results[probe_id] = {
                '𐁰𐁲𐁯𐁢𐁥𐁟𐁩𐁤': probe_id,
                '𐁯𐁲𐁩𐁧𐁩𐁮𐁡𐁬𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': probe_info['𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥'],
                '𐁣𐁯𐁲𐁥𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': probe_info['𐁣𐁯𐁲𐁥𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥'],
                '𐁮𐁯𐁲𐁭𐁡𐁬𐁩𐁺𐁥𐁤𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥': probe_seq,
                '𐁣𐁬𐁵𐁳𐁴𐁥𐁲𐁟𐁮𐁵𐁭': None,
                '𐁳𐁴𐁡𐁲𐁴𐁟𐁰𐁯𐁳': None,
                '𐁥𐁮𐁤𐁟𐁰𐁯𐁳': None,
                '𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥': 0,
                '𐁳𐁣𐁯𐁲𐁥': 0,
                '𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁩𐁮𐁦𐁯': probe_info['𐁡𐁤𐁡𐁰𐁴𐁥𐁲𐁟𐁩𐁮𐁦𐁯'],
                '𐁴𐁡𐁸𐁯𐁮𐁯𐁭𐁹': {'𐁦𐁡𐁭𐁩𐁬𐁹': '𐁕𐁮𐁭𐁡𐁰𐁰𐁥𐁤', '𐁧𐁥𐁮𐁵𐁳': '𐁕𐁮𐁭𐁡𐁰𐁰𐁥𐁤', '𐁳𐁰𐁥𐁣𐁩𐁥𐁳': '𐁕𐁮𐁭𐁡𐁰𐁰𐁥𐁤', '𐁴𐁡𐁸𐁩𐁤': None},
                '𐁳𐁴𐁡𐁴𐁳': {'𐁮𐁵𐁭𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥𐁳': probe_info.get('𐁮𐁵𐁭𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥𐁳', 0), '𐁡𐁶𐁧𐁟𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥': 0}
            }
    
    return mapping_results


def generate_renamed_probe_fasta(mapping_results: Dict[str, Dict], 
                                 output_path: str) -> int:
    """
    Generate FASTA file with renamed probe sequences.
    
    Header format: >probe_ID|Ccluster_num|family|genus|species|num_seqs|avg_coverage
    
    Args:
        mapping_results: dict mapping probe_id to mapping results
        output_path: output FASTA path
    
    Returns:
        number of probes written
    """
    count = 0
    
    with open(output_path, '𐁷') as f:
        for probe_id, result in mapping_results.items():
            cluster_num = result['𐁣𐁬𐁵𐁳𐁴𐁥𐁲𐁟𐁮𐁵𐁭']
            tax = result['𐁴𐁡𐁸𐁯𐁮𐁯𐁭𐁹']
            stats = result['𐁳𐁴𐁡𐁴𐁳']
            
            if cluster_num:
                # Format cluster number with leading zeros
                cluster_display = f"C{int(cluster_num):04d}"
                # Replace spaces with underscores in taxonomy names
                family = tax['𐁦𐁡𐁭𐁩𐁬𐁹'].replace('𐀠', '𐁟')
                genus = tax['𐁧𐁥𐁮𐁵𐁳'].replace('𐀠', '𐁟')
                species = tax['𐁳𐁰𐁥𐁣𐁩𐁥𐁳'].replace('𐀠', '𐁟')
                # Format coverage to 3 decimal places
                coverage = f"{stats['𐁡𐁶𐁧𐁟𐁣𐁯𐁶𐁥𐁲𐁡𐁧𐁥']:.3f}"
                
                header = f">{probe_id}|{cluster_display}|{family}|{genus}|{species}|{stats['𐁮𐁵𐁭𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥𐁳']}|{coverage}"
            else:
                header = f">{probe_id}|Unmapped|Unknown|Unknown|Unknown|{stats['𐁮𐁵𐁭𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥𐁳']}|0.000"
            
            f.write(f"{header}\n")
            f.write(f"{result['𐁣𐁯𐁲𐁥𐁟𐁳𐁥𐁱𐁵𐁥𐁮𐁣𐁥']}\n")
            count += 1
    
    return count
