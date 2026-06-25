#!/usr/bin/env python3
"""
Config module for OligoPlex.

Handles configuration loading, saving, and validation for probe design pipeline.
"""

import json
import os
from typing import Dict, Any, Optional


DEFAULT_CONFIG = {
    'adapter_sequence': 'GTGGAGGTCGCTAGATGGTC',
    'iupac_to_n': True,
    'min_coverage_threshold': 0.80,
    'max_gap_threshold': 0.10,
    'entrez_email': 'user@example.com',
    'entrez_api_key': None,
    'batch_size': 100,
    'threads': None,
    'output_format': 'excel',
    'visualization_width': 1200,
    'visualization_height': 600,
    'metadata_path': None,
}


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from JSON file or use defaults.
    
    Args:
        config_path: path to config JSON file (optional)
    
    Returns:
        dict with configuration parameters
    """
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config = json.load(f)
        config = {**DEFAULT_CONFIG, **config}
    else:
        config = DEFAULT_CONFIG.copy()
    
    return config


def save_config(config: Dict[str, Any], config_path: str) -> None:
    """
    Save configuration to JSON file.
    
    Args:
        config: configuration dict
        config_path: output path for config file
    """
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate configuration parameters.
    
    Args:
        config: configuration dict
    
    Returns:
        True if valid, False otherwise
    """
    errors = []
    
    if config['min_coverage_threshold'] < 0 or config['min_coverage_threshold'] > 1:
        errors.append(f"min_coverage_threshold must be between 0 and 1, got {config['min_coverage_threshold']}")
    
    if config['max_gap_threshold'] < 0 or config['max_gap_threshold'] > 1:
        errors.append(f"max_gap_threshold must be between 0 and 1, got {config['max_gap_threshold']}")
    
    if len(config['adapter_sequence']) < 10:
        errors.append(f"adapter_sequence should be at least 10 bases, got {len(config['adapter_sequence'])}")
    
    if errors:
        print("Configuration validation errors:")
        for error in errors:
            print(f"  - {error}")
        return False
    
    return True
