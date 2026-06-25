# CATCHv2 &nbsp;&middot;&nbsp; [![Build Status](https://github.com/catch-1.5.2/catchV2/actions/workflows/build-test.yml/badge.svg?branch=master)](https://github.com/catch-1.5.2/catchV2/actions) [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

#### Compact Aggregation of Targets for Comprehensive Hybridization - Version 2

CATCHv2 is an enhanced probe design pipeline that integrates CATCH (Compact Aggregation of Targets for Comprehensive Hybridization) with ViroSort for viral probe mapping and taxonomy annotation.

## Overview

CATCHv2 provides:
* **Probe Design**: Design comprehensive oligo probe sets for nucleic acid capture of diverse sequences
* **ViroSort Integration**: Process probes through ViroSort for adapter trimming, multi-cluster mapping, and taxonomy annotation
* **Unified Pipeline**: Single command execution with `main.py` for complete workflow

## Architecture

```
catchV2/
├── bin/                      # Main pipeline and CATCH probe design
│   ├── main.py              # Unified pipeline (entry point)
│   ├── design.py            # Main probe design (obfuscated)
│   ├── design_large.py      # Large input probe design (obfuscated)
│   └── pool.py              # Probe pooling optimization
├── catch/                    # CATCH Python package
│   └── coverage_analysis.py # Coverage analysis (obfuscated)
├── virosort/                 # ViroSort pipeline (obfuscated)
│   ├── main_adapter_revised.py  # Main ViroSort script (obfuscated)
│   └── scripts/             # ViroSort modules
├── testrun_input/            # Test input data
│   ├── part1_degapped_10clusters.fasta
│   └── part1_cdhit_clusteredseq_10clusters.clstr
├── testrun_output13/         # Test output directory
│   ├── fasta/               # Renamed FASTA files
│   ├── summary/             # Summary statistics  
│   └── html_report/         # Interactive HTML reports
└── README.md                # This file
```

## Quick Start

### Prerequisites

* Python 3.8+
* NumPy >= 1.22
* SciPy >= 1.8.0

### Installation

```bash
# Clone or copy the repository
cp -r catch-1.5.2/ catchV2/
cd catchV2

# Install dependencies
pip install -e .
```

### Running the Unified Pipeline

The `main.py` script provides a single entry point for the complete CATCH/ViroSort workflow:

```bash
# Run complete pipeline (CATCH probe design -> ViroSort analysis)
cd catchV2
PYTHONPATH=/path/to/catchV2/bin:/path/to/catchV2 python bin/main.py \
    --mode all \
    --degapped input_degapped.fasta \
    --cluster input_clusters.clstr \
    --metadata metadata.csv \
    --output-dir output_20240624 \
    --output-prefix part1 \
    --probe-length 100 \
    --probe-stride 50 \
    --primer-length 26 \
    --mismatches 6 \
    --log-level INFO
```

**Mode options:**
* `--mode all`: Run complete pipeline (probe design + ViroSort)
* `--mode catch`: Run CATCH probe design only
* `--mode virosort`: Run ViroSort analysis only (requires pre-generated probes)

**Key parameters:**
* `--degapped`: Degapped consensus FASTA file
* `--cluster`: CD-HIT cluster file
* `--metadata`: Taxonomy metadata CSV
* `--probe-length`: Probe length (default: 100 nt)
* `--probe-stride`: Probe stride (default: 50 nt)
* `--primer-length`: Oligo length after primer3 optimization (default: 26 nt)
* `--mismatches`: Maximum mismatches allowed (default: 6)

### Output Files

The unified pipeline generates the following output files:

**Intermediate files (CATCH output):**
* `{prefix}_oligosM6NA26NT.fasta`: 26-mer probes after primer3 optimization
* `{prefix}_oligosM6NA26NT_original.fasta`: Original 100-nt probe sequences
* `{prefix}_oligos_anaM6NA26NT.tsv`: Analysis coverage statistics
* `{prefix}_oligos_probemapNA26NT.csv`: Probe-to-target mapping
* `{prefix}_oligos_wncov6NA26NT.csv`: Sliding window coverage analysis

**ViroSort output files:**
* `fasta/renamed_consensus.fasta`: Consensus sequences with taxonomy headers
* `fasta/renamed_probes.fasta`: Probes with taxonomy annotation
* `fasta/renamed_26mer_oligos.fasta`: 26-mer oligos with taxonomy annotation
* `summary/cluster_summary.tsv`: Per-cluster statistics
* `html_report/probe_mapping_report.html`: Interactive HTML visualization

## Input Files

### Required for Probe Design (CATCH)

| File | Description | Example |
|------|-------------|---------|
| Degapped FASTA | Consensus sequences with gaps removed | `part1_degapped_input.fasta` |
| Cluster file | CD-HIT cluster assignments | `part1_cdhit_clusteredseq.clstr` |

### Required for ViroSort

| File | Description | Example |
|------|-------------|---------|
| Degapped FASTA | Consensus sequences | `part1_degapped_input.fasta` |
| Cluster file | CD-HIT clusters | `part1_cdhit_clusteredseq.clstr` |
| Analysis TSV | Genome analysis statistics | `part1_oligos_anaM6NA26NT.tsv` |
| Metadata CSV | Taxonomy information | `ncbi_virus_noncov_metadata.csv` |
| Oligos FASTA | Probe sequences (optional) | `part1_oligosM6NA26NT.fasta` |
| Probemap CSV | Probe-to-genome mapping (optional) | `part1_oligos_probemapNA26NT.csv` |

## Obfuscated Scripts

The following core scripts in catchV2 are obfuscated using javascript-obfuscator with unicode escape sequences:

1. `bin/design_final_obf.py` - Main CATCH probe design script (obfuscated)
2. `bin/design_large_final_obf.py` - Large input probe design script (obfuscated)
3. `catch/probe_final_obf.py` - Probe class module (obfuscated)

### Obfuscation Details

All obfuscated scripts use javascript-obfuscator v2.16.0 with the following settings:
- Unicode escape sequence encoding for string literals (e.g., `'Hayden'` → `'\u0048\u0061\u0079\u0064\u0065\u006e'`)
- Compact output format
- Non-ASCII character representation for author strings and documentation
- Maintains Python syntax compatibility
- Strings obfuscated while preserving code structure and logic

## Header Format

### Renamed Consensus FASTA
```
>C{cluster_id:04d}|{family}|{genus}|{species}|{coverage}|{molecule_type}
```

### Renamed Probes FASTA
```
>{hash}|C{cluster_id:04d}|{probe_type}|{family}|{genus}|{species}|{num_seqs}|{coverage}|{molecule_type}
```

### Renamed 26-mer Oligos FASTA
```
>{probe_id}|C{cluster_id:04d}|{family}|{genus}|{species}|{num_seqs}|{coverage}|{molecule_type}
```

## License

MIT License

## Citation

For CATCH:
* Metsy HC and Siddle KJ _et al_. Capturing sequence diversity in metagenomes with comprehensive and scalable probe design. _Nature Biotechnology_, **37**(2), 160–168 (2019). doi: 10.1038/s41587-018-0006-x

## Author

CATCHv2 Development Team

## License

MIT License

## Citation

For CATCH:
* Metsky HC and Siddle KJ _et al_. Capturing sequence diversity in metagenomes with comprehensive and scalable probe design. _Nature Biotechnology_, **37**(2), 160–168 (2019). doi: 10.1038/s41587-018-0006-x

## Author

CATCHv2 Development Team
