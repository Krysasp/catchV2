#!/usr/bin/env python3
"""
Taxonomy module for OligoPlex.

Handles taxonomy retrieval using CSV metadata as primary source with NCBI Entrez fallback.
"""

import csv
import os
import time
from typing import Dict, List, Any, Optional


def load_csv_metadata(metadata_path: str) -> Dict[str, Dict[str, str]]:
    """
    Load taxonomy metadata from CSV file indexed by Accession ID.
    
    Args:
        metadata_path: path to CSV file with columns: Accession, Organism_Name, Species, Genus, Family, etc.
    
    Returns:
        dict mapping Accession ID to taxonomy info dict
    """
    metadata_map = {}
    
    if not os.path.exists(metadata_path):
        print(f"  Metadata file not found: {metadata_path}")
        return metadata_map
    
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                accession = row.get('Accession', '').strip()
                if accession:
                    metadata_map[accession] = {
                        'taxid': None,  # Not available in CSV
                        'family': row.get('Family', '').strip() or 'Unknown',
                        'genus': row.get('Genus', '').strip() or 'Unknown',
                        'species': row.get('Species', '').strip() or 'Unknown',
                        'species_display': row.get('Organism_Name', '').strip() or '',
                        'country': row.get('Geo_Location', '').strip() or '',
                        'collection_date': row.get('Collection_Date', '').strip() or '',
                        'molecule_type': row.get('Molecule_type', '').strip() or 'Unknown'
                    }
        print(f"  Loaded metadata for {len(metadata_map)} accessions from CSV")
    except Exception as e:
        print(f"  Error loading metadata CSV: {e}")
    
    return metadata_map


def _get_taxonomy_from_csv(accession: str, metadata_map: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
    """
    Get taxonomy information from pre-loaded CSV metadata.
    
    Args:
        accession: accession ID
        metadata_map: pre-loaded metadata dictionary
    
    Returns:
        dict with taxonomy info, using 'Unknown' for missing fields
    """
    if accession in metadata_map:
        return metadata_map[accession]
    
    # Accession not found - return Unknown defaults
    return {
        'taxid': None,
        'family': 'Unknown',
        'genus': 'Unknown',
        'species': 'Unknown',
        'species_display': '',
        'country': '',
        'collection_date': ''
    }


def _parse_ncbi_response(entry: Any) -> Dict[str, Any]:
    """
    Parse NCBI efetch response into taxonomy dict.
    
    Args:
        entry: parsed XML entry from Entrez.read()
    
    Returns:
        dict with taxid, family, genus, species, country, collection_date
    """
    taxonomy = {
        'taxid': None,
        'family': '',
        'genus': '',
        'species': '',
        'species_display': '',
        'country': '',
        'collection_date': ''
    }
    
    taxonomy['species_display'] = entry.get('GBSeq_definition', '')
    taxonomy['species'] = entry.get('GBSeq_organism', '')
    
    # Extract taxid from db_xref
    for feature in entry.get('GBSeq_feature-table', []):
        if feature.get('GBFeature_key') == 'source':
            for qualifier in feature.get('GBFeature_quals', []):
                qual_name = qualifier.get('GBQualifier_name', '')
                if qual_name == 'db_xref':
                    db_value = qualifier.get('GBQualifier_value', '')
                    if db_value.startswith('taxon:'):
                        taxonomy['taxid'] = int(db_value.split(':')[1])
                elif qual_name == 'country':
                    taxonomy['country'] = qualifier.get('GBQualifier_value', '')
                elif qual_name == 'collection_date':
                    taxonomy['collection_date'] = qualifier.get('GBQualifier_value', '')
    
    # Parse taxonomy lineage from GBSeq_taxonomy
    taxonomy_lineage = entry.get('GBSeq_taxonomy', '')
    if taxonomy_lineage:
        lineage_parts = [part.strip() for part in taxonomy_lineage.split(';')]
        for part in lineage_parts:
            if not taxonomy['family'] and part.endswith('viridae'):
                taxonomy['family'] = part
            elif not taxonomy['genus'] and part.endswith('virus'):
                taxonomy['genus'] = part
    
    # If family/genus not found, try organism name
    if not taxonomy['family'] or not taxonomy['genus']:
        organism = entry.get('GBSeq_organism', '')
        if organism:
            family_patterns = [p for p in taxonomy_lineage.split(';') if p.strip().endswith('viridae')]
            if family_patterns:
                taxonomy['family'] = family_patterns[0].strip()
            
            words = organism.split()
            if len(words) >= 2:
                taxonomy['genus'] = words[0]
                taxonomy['species'] = ' '.join(words[:2])
    
    return taxonomy


def get_taxonomy_from_accession(accession: str, config: Dict[str, Any], metadata_map: Optional[Dict[str, Dict[str, str]]] = None) -> Dict[str, Any]:
    """
    Get taxonomy information for a single accession.
    
    Primary source: CSV metadata
    Fallback: NCBI Entrez if metadata not available
    
    Args:
        accession: accession ID (e.g., 'NC_129209')
        config: configuration dict with entrez settings and metadata path
        metadata_map: optional pre-loaded metadata dictionary
    
    Returns:
        dict with taxid, family, genus, species, country, collection_date
    """
    # Try CSV metadata first
    if metadata_map is not None:
        if accession in metadata_map:
            return metadata_map[accession]
    
    # CSV fallback: load metadata on demand
    metadata_path = config.get('metadata_path')
    if metadata_path and os.path.exists(metadata_path):
        if metadata_map is None:
            metadata_map = load_csv_metadata(metadata_path)
        if accession in metadata_map:
            return metadata_map[accession]
    
    # NCBI Entrez fallback
    return _fetch_from_ncbi(accession, config)


def _fetch_from_ncbi(accession: str, config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Fetch taxonomy from NCBI Entrez as fallback.
    
    Args:
        accession: accession ID
        config: configuration dict with entrez settings
    
    Returns:
        dict with taxonomy info
    """
    from Bio import Entrez
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            Entrez.email = config.get('entrez_email', 'user@example.com')
            if config.get('entrez_api_key'):
                Entrez.api_key = config['entrez_api_key']
            
            handle = Entrez.efetch('nuccore', id=accession, retmode='xml', rettype='gb')
            response = Entrez.read(handle)
            handle.close()
            
            if response:
                return _parse_ncbi_response(response[0])
            
            time.sleep(1.0 * (attempt + 1))
        except Exception as e:
            error_msg = str(e)
            if 'HTTP Error 429' in error_msg:
                time.sleep(5.0 * (attempt + 1))
            elif 'HTTP Error 400' in error_msg:
                time.sleep(2.0 * (attempt + 1))
            elif attempt < max_retries - 1:
                time.sleep(2.0 * (2 ** attempt))
            else:
                print(f"    NCBI fallback failed for {accession}: {e}")
                return {
                    'taxid': None,
                    'family': 'Unknown',
                    'genus': 'Unknown',
                    'species': 'Unknown',
                    'species_display': '',
                    'country': '',
                    'collection_date': ''
                }
    
    return {
        'taxid': None,
        'family': 'Unknown',
        'genus': 'Unknown',
        'species': 'Unknown',
        'species_display': '',
        'country': '',
        'collection_date': ''
    }


def get_taxonomy_bulk(accessions: List[str], config: Dict[str, Any], metadata_map: Optional[Dict[str, Dict[str, str]]] = None) -> Dict[str, Dict[str, str]]:
    """
    Get taxonomy information for multiple accessions.
    
    Uses CSV metadata as primary source (instant lookup).
    Only accessions not in metadata fall back to NCBI Entrez.
    
    Args:
        accessions: list of accession IDs
        config: configuration dict with entrez settings
        metadata_map: optional pre-loaded metadata dictionary
    
    Returns:
        dict mapping accession to taxonomy info dict
    """
    # Load metadata if not provided
    if metadata_map is None:
        metadata_path = config.get('metadata_path')
        if metadata_path:
            metadata_map = load_csv_metadata(metadata_path)
        else:
            metadata_map = {}
    
    taxonomy_map = {}
    ncbi_accessions = []
    
    # Separate accessions: CSV vs NCBI
    for acc in accessions:
        if acc in metadata_map:
            taxonomy_map[acc] = metadata_map[acc]
        else:
            ncbi_accessions.append(acc)
    
    print(f"  CSV lookup: {len([a for a in accessions if a in metadata_map])} accessions, NCBI fallback: {len(ncbi_accessions)} accessions")
    
    # Fetch missing accessions from NCBI
    if ncbi_accessions:
        print(f"  Fetching {len(ncbi_accessions)} accessions from NCBI...")
        for acc in ncbi_accessions:
            taxonomy_map[acc] = _fetch_from_ncbi(acc, config)
    
    return taxonomy_map


def get_taxonomy_from_cluster_file(cluster_path: str, config: Dict[str, Any], cluster_sequences: Dict[str, List[Dict[str, Any]]] = None) -> Dict[str, Dict[str, str]]:
    """
    Get taxonomy information for all sequences in a cluster file.
    
    Args:
        cluster_path: path to cluster file
        config: configuration dict with entrez settings and metadata_path
        cluster_sequences: optional pre-parsed cluster sequences
    
    Returns:
        dict mapping cluster_num to dict mapping accession to taxonomy info
    """
    from scripts.parsers import parse_cluster_file
    
    # Load metadata once
    metadata_map = None
    metadata_path = config.get('metadata_path')
    if metadata_path:
        metadata_map = load_csv_metadata(metadata_path)
    
    if cluster_sequences is None:
        print("Parsing cluster file...")
        cluster_sequences = parse_cluster_file(cluster_path)
    else:
        print(f"Using pre-parsed cluster file with {len(cluster_sequences)} clusters...")
    
    # Collect all unique accessions
    all_accessions = set()
    for cluster_seqs in cluster_sequences.values():
        for seq_info in cluster_seqs:
            all_accessions.add(seq_info['accession'])
    
    print(f"Retrieving taxonomy for {len(all_accessions)} unique accessions...")
    taxonomy_map = get_taxonomy_bulk(list(all_accessions), config, metadata_map)
    
    # Organize by cluster
    cluster_taxonomy = {}
    for cluster_num, cluster_seqs in cluster_sequences.items():
        cluster_taxonomy[cluster_num] = {}
        for seq_info in cluster_seqs:
            acc = seq_info['accession']
            cluster_taxonomy[cluster_num][acc] = taxonomy_map.get(acc, {
                'taxid': None,
                'family': 'Unknown',
                'genus': 'Unknown',
                'species': 'Unknown',
                'species_display': '',
                'country': '',
                'collection_date': ''
            })
    
    return cluster_taxonomy


def derive_species_name(family: str, genus: str, species_display: str) -> str:
    """
    Derive a standardized species name from taxonomy information.
    
    Args:
        family: family name
        genus: genus name
        species_display: species display name
    
    Returns:
        standardized species name with spaces replaced by underscores
    """
    if species_display:
        species = species_display.replace(' ', '_')
    elif genus and genus != 'Unknown':
        species = genus.replace(' ', '_')
    else:
        species = 'Unknown'
    
    return species


def get_cluster_taxonomy_summary(cluster_taxonomy: Dict[str, Dict[str, Dict]]) -> Dict[str, Dict[str, str]]:
    """
    Summarize taxonomy information by cluster.
    
    Args:
        cluster_taxonomy: dict mapping cluster_num to accession -> taxonomy info
    
    Returns:
        dict mapping cluster_num to cluster-level taxonomy summary
    """
    cluster_summary = {}
    
    for cluster_num, accessions_tax in cluster_taxonomy.items():
        families = set()
        genera = set()
        species_names = set()
        
        for acc_info in accessions_tax.values():
            if acc_info.get('family') and acc_info['family'] != 'Unknown':
                families.add(acc_info['family'])
            if acc_info.get('genus') and acc_info['genus'] != 'Unknown':
                genera.add(acc_info['genus'])
            if acc_info.get('species') and acc_info['species'] != 'Unknown':
                species_names.add(acc_info['species'])
        
        cluster_summary[cluster_num] = {
            'family': list(families)[0] if len(families) == 1 else ';'.join(sorted(families)) if families else 'Unknown',
            'genus': list(genera)[0] if len(genera) == 1 else ';'.join(sorted(genera)) if genera else 'Unknown',
            'species': list(species_names)[0] if len(species_names) == 1 else ';'.join(sorted(species_names)) if species_names else 'Unknown',
            'taxid': list(accessions_tax.values())[0].get('taxid') if accessions_tax else None
        }
    
    return cluster_summary
