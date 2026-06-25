#!/usr/bin/env python3
"""
Parsers module for OligoPlex.

Handles FASTA, cluster, probemap, and analysis file parsing.
"""

import re
from typing import Dict, List, Tuple, Any, Optional


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse FASTA file and return dictionary of sequences.
    
    Args:
        fasta_path: path to FASTA file
    
    Returns:
        dict mapping headers to sequences
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


def parse_cluster_file(cluster_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Parse cluster file and return dictionary mapping cluster numbers to sequence info.
    
    Cluster file format (CD-HIT):
    >Cluster 11
    0   35927aa, >PX778758|Human... at 99.16%
    1   35947aa, >PX778759|Human... at 99.04%
    2   35980aa, >PV125340|Human... at 99.43%
    
    Args:
        cluster_path: path to cluster file
    
    Returns:
        dict mapping cluster_num to list of sequence info dicts with keys:
        - accession: accession ID
        - length: sequence length
        - header: full FASTA header
        - index: position in cluster
    """
    clusters = {}
    current_cluster = None
    
    with open(cluster_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>Cluster '):
                cluster_num = line.split()[1]
                current_cluster = cluster_num
                clusters[cluster_num] = []
            else:
                # Parse cluster member line
                # Format: "0   35927aa, >PX778758|Human... at 99.16%"
                match = re.search(r'>([A-Za-z0-9_.]+)', line)
                if match:
                    accession = match.group(1)
                    
                    length_match = re.search(r'(\d+)aa', line)
                    length = int(length_match.group(1)) if length_match else 0
                    
                    # Extract full header - everything after > until whitespace or ...
                    header_match = re.search(r'>([^\s]+)', line)
                    header = header_match.group(1) if header_match else accession
                    
                    clusters[current_cluster].append({
                        'accession': accession,
                        'length': length,
                        'header': header,
                        'index': len(clusters[current_cluster])
                    })
    
    return clusters


def parse_probemap(probemap_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Parse probemap TSV file containing probe sequences.
    
    Format:
    Probe identifier	Probe sequence	Number sequences mapped to
    6498660b6e	AGGCCCTGGCTGCTGATATGCATCAATAATATACCTTACACTGGATTTGAGCCAATATTAAAATGAAGTGGGCGGAGTGAATAGTTAATTGACCTTTTGGGACAGCGGTG	18
    
    Args:
        probemap_path: path to probemap TSV file
    
    Returns:
        dict mapping probe identifiers to probe info dicts with keys:
        - identifier: probe ID
        - sequence: probe sequence
        - num_sequences: number of sequences mapped to
    """
    probes = {}
    
    with open(probemap_path, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                identifier = parts[0]
                sequence = parts[1].upper()
                num_sequences = int(parts[2])
                
                probes[identifier] = {
                    'identifier': identifier,
                    'sequence': sequence,
                    'num_sequences': num_sequences
                }
    
    return probes


def parse_analysis_file(analysis_path: str) -> Dict[str, Dict[str, float]]:
    """
    Parse analysis TSV file containing coverage statistics.
    
    Format:
    Genome	Num bases covered	Frac bases covered	Frac bases covered over unambig	Average coverage/depth	Average coverage/depth over unambig
    final_consensus_sequences.fasta, genome 0	36511	1.0	1.0	3.6395880693489633	3.6395880693489633
    
    Args:
        analysis_path: path to analysis TSV file
    
    Returns:
        dict mapping cluster numbers to stats dicts with keys:
        - num_bases_covered: number of bases covered
        - frac_bases_covered: fraction of bases covered
        - frac_bases_over_unambig: fraction over unambiguous
        - avg_coverage_depth: average coverage depth
        - avg_coverage_over_unambig: average coverage over unambiguous
    """
    stats = {}
    
    with open(analysis_path, 'r') as f:
        header = f.readline().strip().split('\t')
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) >= 6:
                # Extract cluster number from genome header
                # Format: "final_consensus_sequences.fasta, genome 0"
                genome_header = parts[0]
                cluster_match = re.search(r'genome\s+(\d+)', genome_header)
                
                if cluster_match:
                    cluster_num = cluster_match.group(1)
                    
                    stats[cluster_num] = {
                        'num_bases_covered': int(parts[1]),
                        'frac_bases_covered': float(parts[2]),
                        'frac_bases_over_unambig': float(parts[3]),
                        'avg_coverage_depth': float(parts[4]),
                        'avg_coverage_over_unambig': float(parts[5])
                    }
    
    return stats


def extract_consensus_clusters(consensus_fasta: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    """
    Extract consensus sequences organized by cluster number.
    
    Consensus FASTA header format:
    >CONS|0011|95pct85pct_degap
    
    Where '11' is the cluster number.
    
    Args:
        consensus_fasta: dict mapping FASTA headers to sequences
    
    Returns:
        dict mapping cluster_num to dict with 'header' and 'sequence' keys
    """
    clusters = {}
    
    for header, sequence in consensus_fasta.items():
        # Parse header: CONS|0011|95pct85pct_degap
        # Header format is like "CONS|0011|95pct85pct_degap"
        if '|' in header:
            parts = header.split('|')
            if len(parts) >= 2:
                # Extract cluster number from second part
                cluster_part = parts[1]
                # Remove any leading zeros for comparison
                cluster_num = str(int(cluster_part))
                clusters[cluster_num] = {
                    'header': header,
                    'sequence': sequence
                }
    
    return clusters


def extract_raw_sequences_with_taxonomy(raw_fasta: Dict[str, str], cluster_sequences: Dict[str, List[Dict]]) -> Dict[str, Dict[str, Any]]:
    """
    Extract raw sequences with full metadata from FASTA based on cluster file.
    
    Raw FASTA header format:
    >PX778758|Human_adenovirus_4_|Adenoviridae|Mastadenovirus|Finland|2014-02-24
    
    Args:
        raw_fasta: dict mapping FASTA headers to sequences
        cluster_sequences: dict mapping cluster_num to list of sequence info
    
    Returns:
        dict mapping accession to sequence info with full metadata
    """
    sequences = {}
    
    # Build set of accessions from cluster file
    cluster_accessions = set()
    for cluster_seqs in cluster_sequences.values():
        for seq_info in cluster_seqs:
            cluster_accessions.add(seq_info['accession'])
    
    # Extract matching sequences from raw FASTA
    for header, sequence in raw_fasta.items():
        accession = header.split('|')[0]
        
        if accession in cluster_accessions:
            # Parse full header for metadata
            parts = header.split('|')
            
            metadata = {
                'accession': accession,
                'header': header,
                'sequence': sequence,
                'species': parts[1] if len(parts) > 1 else '',
                'family': parts[2] if len(parts) > 2 else '',
                'genus': parts[3] if len(parts) > 3 else '',
                'location': parts[4] if len(parts) > 4 else '',
                'date': parts[5] if len(parts) > 5 else ''
            }
            
            sequences[accession] = metadata
    
    return sequences


def extract_raw_sequences(raw_fasta: Dict[str, str], cluster_sequences: Dict[str, List[Dict]]) -> Dict[str, Dict[str, Any]]:
    """
    Extract raw sequences with full metadata from FASTA based on cluster file.
    
    Raw FASTA header format:
    >PX778758|Human_adenovirus_4_|Adenoviridae|Mastadenovirus|Finland|2014-02-24
    
    Args:
        raw_fasta: dict mapping FASTA headers to sequences
        cluster_sequences: dict mapping cluster_num to list of sequence info
    
    Returns:
        dict mapping accession to sequence info with full metadata
    """
    sequences = {}
    
    # Build set of accessions from cluster file
    cluster_accessions = set()
    for cluster_seqs in cluster_sequences.values():
        for seq_info in cluster_seqs:
            cluster_accessions.add(seq_info['accession'])
    
    # Extract matching sequences from raw FASTA
    for header, sequence in raw_fasta.items():
        accession = header.split('|')[0]
        
        if accession in cluster_accessions:
            # Parse full header for metadata
            parts = header.split('|')
            
            metadata = {
                'accession': accession,
                'header': header,
                'sequence': sequence,
                'species': parts[1] if len(parts) > 1 else '',
                'family': parts[2] if len(parts) > 2 else '',
                'genus': parts[3] if len(parts) > 3 else '',
                'location': parts[4] if len(parts) > 4 else '',
                'date': parts[5] if len(parts) > 5 else ''
            }
            
            sequences[accession] = metadata
    
    return sequences
