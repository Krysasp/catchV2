# CATCHv2 &nbsp;&middot;&nbsp; [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

#### Compact Aggregation of Targets for Comprehensive Hybridization - Version 2

CATCHv2 is an enhanced probe design pipeline that integrates CATCH (Compact Aggregation of Targets for Comprehensive Hybridization) with ViroSort for viral probe mapping and taxonomy annotation.

## Usage, Utility, and Functionalities

CATCHv2 is a comprehensive bioinformatics pipeline designed for viral probe development and analysis. It combines CATCH for probe design and ViroSort as helper tool, created in-house, for viral sequence characterization—into a unified workflow.

## Examples of oPools 
Find the following links for oPools designed from complete genome sequences from Genbank
1. Human-related respiratory viruses (non-EV) (https://thinfi.com/0lvhy)
2. Other human viruses (https://thinfi.com/0lvi0)
3. Human and other mammals (zoonotic or broad mammalian host range) (https://thinfi.com/0lvi1)
4. Primarily non-human vertebrates (https://thinfi.com/0lvi2)
5. Avian viruses (transmissible to humans) (https://thinfi.com/0lvi3)
6. Plant-associated viruses (https://thinfi.com/0lvi4)
7. Fungi-associated (mycoviruses) (https://thinfi.com/0lvi5)
8. Bacteriophages (https://thinfi.com/0lvi7)
9. Environmental viruses (commonly detected in wastewater, sewage, soil, feces, or environmental samples; host often unknown or likely fungi/protists/invertebrates) (https://thinfi.com/0lvi8)


### Primary Use Cases

**1. Viral Probe Design for NGS Enrichment & Target Amplification**
- Design oligonucleotide probes (baits) that hybridize to diverse viral sequences
- Capture sequence diversity from metagenomic samples or viral populations
- Optimize probes for maximum coverage across viral families, genera, and species
- Enable reselection of sequence frames from probes to serve as spike primers 
- Handle highly variable viral genomes with degenerate probe design

**2. Multi-Cluster Probe Pooling**
- Aggregate probes from multiple viral clusters into optimized pools
- Balance probe representation to ensure even hybridization across targets
- Minimize redundancy while maintaining comprehensive coverage
- Design primers (iether in forward/reverse orientation) flanking probe regions

**3. Viral Taxonomy Annotation**
- Annotate probe sequences with taxonomic information (family, genus, species)
- Map probes to reference genomes and track coverage statistics
- Generate taxonomy-aware FASTA headers for downstream analysis
- Identify probe specificity and cross-reactivity patterns

**4. Quality Control and Coverage Analysis**
- Calculate sliding window coverage across reference sequences
- Identify poorly covered regions (gaps) in probe design
- Assess probe-to-genome mapping quality
- Generate comprehensive statistics for probe set evaluation

### How It Works (Workflow)

```
Input (viral sequences)
    │
    ├── Step 1: Clustering (CD-HIT)
    │       └── Group similar sequences into clusters
    │
    ├── Step 2: Consensus Generation
    │       └── Create degapped consensus per cluster
    │
    ├── Step 3: Probe Design (CATCH)
    │       ├── Extract #1-nt probe sequences at specified stride
    │       ├── Optimize #2-mer oligos using primer3
    │       └── Filter probes by mismatch tolerance
    │
    ├── Step 4: Probe Mapping (ViroSort)
    │       ├── Map probes to reference genomes
    │       ├── Annotate with taxonomy metadata
    │       └── Calculate coverage statistics
    │
    └── Output (annotated probe sets + reports)
```

### Key Functionalities

| Function | Description | Output Impact |
|----------|-------------|---------------|
| **Degapped Consensus** | Removes gaps from aligned sequences to create continuous probe templates | Ensures probes target conserved regions |
| **CD-HIT Integration** | Clusters sequences by similarity threshold (default: 90%) | Reduces redundancy, groups related viruses |
| **Probe Extraction** | Sliding window extraction of probe sequences | Generates comprehensive probe library |
| **Primer3 Optimization** | Designs flanking primers for PCR amplification | Enables probe synthesis via PCR |
| **Mismatch Tolerance** | Allows up to N mismatches per probe | Accommodates viral sequence diversity |
| **Taxonomy Annotation** | Adds family/genus/species to FASTA headers | Enables downstream classification |
| **Coverage Analysis** | Calculates per-base and per-region coverage | Identifies gaps in probe coverage |
| **HTML Visualization** | Generates interactive probe mapping reports | Facilitates manual inspection |

## Architecture

The CATCHv2 project is organized into modular components that handle distinct stages of the probe design and analysis workflow:

```
catchV2/
├── bin/                      # Main pipeline and CATCH probe design
│   ├── main.py              # Unified pipeline orchestrator
│   │                        # • Coordinates CATCH and ViroSort execution
│   │                        # • Manages input validation and file I/O
│   │                        # • Handles error recovery and logging
│   ├── design.py            # CATCH core probe design engine (obfuscated)
│   │                        # • Extracts probe sequences from consensus
│   │                        # • Calculates probe diversity metrics
│   │                        # • Applies mismatch tolerance filtering
│   ├── design_large.py      # Optimized for large input datasets (obfuscated)
│   │                        # • Memory-efficient processing of thousands of clusters
│   │                        # • Batch processing with chunked I/O
│   └── pool.py              # Probe pooling and primer design
│                            # • Designs flanking primers (forward/reverse)
│                            # • Optimizes probe pool composition
│                            # • Balances probe representation
│
├── catch/                    # CATCH Python package
│   ├── probe.py             # Probe class and sequence manipulation (obfuscated)
│   │                        # • Represents probe objects with metadata
│   │                        # • Implements sequence operations (reverse complement, etc.)
│   │                        # • Handles FASTA parsing and output
│   ├── coverage_analysis.py # Coverage statistics and gap detection (obfuscated)
│   │                        # • Calculates sliding window coverage
│   │                        # • Identifies uncovered regions
│   │                        # • Generates coverage heatmaps
│   └── utils/               # Utility modules
│       ├── seq_io.py        # FASTA/FASTQ parsing and writing
│       ├── interval.py      # Genomic interval operations
│       └── longest_common_substring.py
│
├── virosort/                 # ViroSort viral analysis pipeline (obfuscated)
│   ├── main_adapter_revised.py  # Main ViroSort orchestrator
│   │                        # • Adapter trimming of probe sequences
│   │                        # • Multi-cluster probe mapping
│   │                        # • Taxonomy annotation pipeline
│   └── scripts/             # ViroSort modules
│       ├── config.py        # Configuration and parameter management
│       ├── parsers.py       # FASTA/TSV/CSV parsing utilities
│       ├── probe_mapper.py  # Probe-to-genome alignment
│       ├── taxonomy.py      # Taxonomy lookup and annotation
│       ├── output_writer.py # Formatted output generation
│       └── visualizer.py    # HTML report generation
│
├── testrun_input/            # Test input data (example degapped FASTA and CD-HIT cluster files)
│
├── testrun_output13/         # Example output directory
│   ├── fasta/               # Renamed FASTA files (consensus, probes, oligos)
│   ├── summary/             # Summary statistics (cluster_summary.tsv)
│   └── html_report/         # Interactive HTML reports
│
└── README.md                # This documentation file
```

### Component Interactions

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   bin/main.py   │────▶│  bin/design.py  │────▶│   catch/probe.py│
│   (orchestrator)│     │ (probe design)  │     │ (probe objects) │
└────────┬────────┘     └────────┬────────┘     └────────┬────────┘
         │                       │                        │
         │         ┌─────────────┴─────────────┐        │
         │         │   catch/coverage_analysis │        │
         │         │   (coverage statistics)   │        │
         │         └─────────────┬─────────────┘        │
         │                       │                      │
         ▼                       ▼                      ▼
┌─────────────────────────────────────────────────────────────┐
│                    bin/virosort/main.py                     │
│         (adapter trimming, mapping, taxonomy annotation)    │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
         ┌────────────────────┴────────────────────┐
         ▼                                         ▼
┌─────────────────────┐              ┌─────────────────────────┐
│  Output FASTA files │              │  Output TSV/HTML reports│
│  (annotated probes) │              │  (statistics, mapping)  │
└─────────────────────┘              └─────────────────────────┘
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

The unified pipeline generates comprehensive output files organized into three categories:

#### 1. Intermediate Files (CATCH Output)

| File | Format | Description | Usage |
|------|--------|-------------|-------|
| `{prefix}_oligomers.fasta` | FASTA | 26-mer oligo probes after primer3 optimization | Synthesis-ready probe sequences |
| `{prefix}_oligos_original.fasta` | FASTA | Original 100-nt probe sequences before optimization | Reference for probe design validation |
| `{prefix}_oligos_ana.tsv` | TSV | Analysis coverage statistics | Quality control, coverage assessment |
| `{prefix}_oligos_probemap.csv` | CSV | Probe-to-target genome mapping | Specificity analysis, off-target detection |
| `{prefix}_oligos_wncov.csv` | CSV | Sliding window coverage analysis | Gap identification, coverage uniformity |

#### 2. Example of Output Files

**FASTA Directory (`fasta/`):**

| File | Description | Header Format |
|------|-------------|---------------|
| `renamed_consensus.fasta` | Consensus sequences with taxonomy annotation | `>C{cluster_id}|{family}\|{genus}\|{species}\|{coverage}\|{molecule_type}` |
| `renamed_probes.fasta` | 100-nt probes with full taxonomy metadata | `>{hash}|C{cluster_id}\|{probe_type}\|{family}\|{genus}\|{species}\|{num_seqs}\|{coverage}\|{molecule_type}` |
| `renamed_26mer_oligos.fasta` | Optimized 26-mer oligos with annotation | `>{probe_id}|C{cluster_id}\|{family}\|{genus}\|{species}\|{num_seqs}\|{coverage}\|{molecule_type}` |

**Summary Directory (`summary/`):**

| File | Description | Key Metrics |
|------|-------------|-------------|
| `cluster_summary.tsv` | Per-cluster statistics | Cluster ID, sequence count, probe count, coverage percentage, taxonomy hierarchy |

**HTML Report Directory (`html_report/`):**

| File | Description | Interactive Features |
|------|-------------|---------------------|
| `probe_mapping_report.html` | Visualization of probe-to-genome mapping | Clickable probe entries, coverage heatmaps, taxonomy filters |

#### 3. Output File Summaries

**Coverage Statistics (`*_ana*.tsv`):**
```
genome_id    total_probes    covered_bases    coverage_pct    avg_mismatches
C0001        45              4420             95.3            1.2
```
- `total_probes`: Number of probes designed for this genome
- `covered_bases`: Unique bases covered by at least one probe
- `coverage_pct`: Percentage of genome covered
- `avg_mismatches`: Average mismatches per probe-genome alignment

**Probe Mapping (`*_probemap*.csv`):**
```
probe_id    target_genome    start    end    mismatches    probe_type
probe_001   C0001            156      255    3             cluster_specific
```
- `probe_id`: Unique probe identifier (hash-based)
- `target_genome`: Reference genome cluster ID
- `start/end`: Genomic coordinates of probe alignment
- `mismatches`: Number of mismatches in alignment
- `probe_type`: Classification (cluster_specific, cross_cluster, etc.)

**Window Coverage (`*_wncov*.csv`):**
```
genome_id    window_start    window_end    probes_in_window    coverage_status
C0001        1               100           2                   covered
C0001        101             200           0                   uncovered
```
- `window_start/end`: Genomic window coordinates
- `probes_in_window`: Number of probes overlapping this window
- `coverage_status`: `covered` (≥1 probe) or `uncovered` (0 probes)

## Input Files

### Required for Probe Design (CATCH)

| File | Description | Example |
|------|-------------|---------|
| Degapped FASTA | Consensus sequences with gaps removed | `part1_degapped_input.fasta` |
| Cluster file | CD-HIT cluster assignments | `part1_cdhit_clusteredseq.clstr` |

### Required for ViroSort

| File | Description | Example |
|------|-------------|---------|
| Degapped FASTA | Consensus sequences | `*_input.fasta` |
| Cluster file | CD-HIT clusters | `*_clusteredseq.clstr` |
| Analysis TSV | Genome analysis statistics | `*_analysis.tsv` |
| Metadata CSV | Taxonomy information | `*_metadata.csv` |
| Oligos FASTA | Probe sequences (optional) | `*_probes.fasta` |
| Probemap CSV | Probe-to-genome mapping (optional) | `*_probemap.csv` |

### Obfuscation Details
The following scripts have been obfuscated:
- virosort/main_adapter_revised.py
- virosort/probe_mapper.py
- virosort/taxonomy.py
- bin/design_large.py
- bin/design.py
  
### Email to 'jonatasp92@gmail.com' for more details.

## License

MIT License

## Citation
Initial design of CATCH:
Metsky, H. C., Siddle, K. J., et. al. (2019). Capturing sequence diversity in metagenomes with comprehensive and scalable probe design. Nature biotechnology, 37(2), 160–168. https://doi.org/10.1038/s41587-018-0006-x

## Author
Jonathan Chan, University of Malaysia Sarawak
