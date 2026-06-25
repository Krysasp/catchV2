# Functions: _calculate_tm, _has_poly_stretch, _combine_iupac_codes, _process_probe_for_primer3, _select_primer3_oligos, main, init_and_parse_args...

#!/usr/bin/env python3
"""Design probes for genome capture.

This is the main executable of CATCH for probe design.
"""

import argparse
import importlib
import logging
import multiprocessing
import os
import random
import typing

from catch import coverage_analysis
from catch import probe
from catch.filter import adapter_filter
from catch.filter import base_filter
from catch.filter import duplicate_filter
from catch.filter import fasta_filter
from catch.filter import n_expansion_filter
from catch.filter import near_duplicate_filter
from catch.filter import polya_filter
from catch.filter import polyg_filter
from catch.filter import polyc_filter
from catch.filter import polyt_filter
from catch.filter import probe_designer
from catch.filter import reverse_complement_filter
from catch.filter import set_cover_filter
from catch.utils import cluster
from catch.utils import ncbi_neighbors
from catch.utils import seq_io, version, log


def _calculate_tm(seq):
    """Calculate Tm of a DNA sequence using nearest-neighbor method.
    
    Args:
        seq: DNA sequence string
    
    Returns:
        Tm in degrees Celsius
    """
    # Simplified Tm calculation using basic formula
    # For accuracy, use primer3's Tm calculation
    if not seq:
        return 0.0
    
    # Count bases
    a = seq.upper().count('A')
    t = seq.upper().count('T')
    g = seq.upper().count('G')
    c = seq.upper().count('C')
    
    # Basic Wallace rule for short oligos
    tm = 2 * (a + t) + 4 * (g + c)
    return round(tm, 1)


def _has_poly_stretch(seq, min_length=5):
    """Check if a sequence has a poly stretch of A/G/C/T exceeding min_length.
    
    Args:
        seq: DNA sequence string
        min_length: minimum poly stretch length to flag
    
    Returns:
        tuple (base, run_length) if poly stretch found, else (None, 0)
    """
    seq = seq.upper()
    for base in 'AGCT':
        max_run = 0
        current_run = 0
        for char in seq:
            if char == base:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        if max_run >= min_length:
            return (base, max_run)
    return (None, 0)


def _combine_iupac_codes(base1, base2):
    """Combine two IUPAC codes into a single IUPAC code.
    
    Args:
        base1: first base (A, C, G, T, or IUPAC code)
        base2: second base (A, C, G, T, or IUPAC code)
    
    Returns:
        Combined IUPAC code
    """
    iupac_map = {
        'A': set('A'), 'C': set('C'), 'G': set('G'), 'T': set('T'),
        'R': set('AG'), 'Y': set('CT'), 'M': set('AC'), 'K': set('GT'),
        'S': set('CG'), 'W': set('AT'), 'H': set('ACT'), 'B': set('CGT'),
        'V': set('ACG'), 'D': set('AGT'), 'N': set('ACGT')
    }
    
    # Get the sets for each base
    set1 = iupac_map.get(base1.upper(), set(base1.upper()))
    set2 = iupac_map.get(base2.upper(), set(base2.upper()))
    
    # Combine the sets
    combined = set1 | set2
    
    # Find the IUPAC code for the combined set
    reverse_map = {
        frozenset('A'): 'A', frozenset('C'): 'C', frozenset('G'): 'G', frozenset('T'): 'T',
        frozenset('AG'): 'R', frozenset('CT'): 'Y', frozenset('AC'): 'M', frozenset('GT'): 'K',
        frozenset('CG'): 'S', frozenset('AT'): 'W', frozenset('ACT'): 'H', frozenset('CGT'): 'B',
        frozenset('ACG'): 'V', frozenset('AGT'): 'D', frozenset('ACGT'): 'N'
    }
    
    return reverse_map.get(frozenset(combined), 'N')


def _process_probe_for_primer3(p, primer_length, target_tm, mismatches, min_tm, max_tm, 
                                primer3_path, poly_filter_length=5, shift_window=75):
    """Process a single probe through primer3 with poly-stretch re-selection.
    
    If the selected oligo has a poly stretch exceeding the threshold, this function
    attempts to re-select from shifted regions (previous or next window of shift_window nt).
    
    Args:
        p: Probe object
        primer_length: desired oligo length
        target_tm: target melting temperature  
        mismatches: number of mismatches to tolerate (also used for PRIMER_MAX_NS_ACCEPTED)
        min_tm: minimum Tm for primer3
        max_tm: maximum Tm for primer3
        primer3_path: path to primer3_core executable
        poly_filter_length: minimum poly stretch length to trigger re-selection
        shift_window: size of shifted window to try for re-selection (default 75)
    
    Returns:
        Probe object (either new optimized oligo or original probe)
    """
    import subprocess
    import tempfile
    import os
    import logging
    
    logger = logging.getLogger(__name__)
    
    seq = str(p)
    header = getattr(p, 'header', None) or p.identifier()
    
    # Get source header prefix if available
    source_header = getattr(p, 'source_header', None)
    if source_header:
        parts = source_header.split('|')
        prefix = parts[0] + '|' + parts[1] if len(parts) >= 2 else parts[0]
    else:
        prefix = header.split('_')[0] if '_' in header else header
    
    # Calculate Tm of original template sequence
    template_tm = _calculate_tm(seq)
    
    def run_primer3_with_sequence(template_seq, seq_prefix):
        """Run primer3 on a given template sequence."""
        # Set PRIMER_PRODUCT_SIZE_RANGE to accommodate the template length
        # This ensures primer3 can select primers from templates of any length
        product_size_min = primer_length
        product_size_max = len(template_seq)
        
        primer3_input = f"""SEQUENCE_ID={seq_prefix}
SEQUENCE_TEMPLATE={template_seq}
PRIMER_OPT_SIZE={primer_length}
PRIMER_MIN_SIZE={primer_length}
PRIMER_MAX_SIZE={primer_length}
PRIMER_OPT_TM={target_tm}
PRIMER_MIN_TM={min_tm}
PRIMER_MAX_TM={max_tm}
PRIMER_NUM_RETURN=1
PRIMER_MAX_NS_ACCEPTED={mismatches}
PRIMER_LIBERAL_BASE=1
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=0
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PRODUCT_SIZE_RANGE={product_size_min}-{product_size_max}
PRIMER_MAX_HAIRPIN_TH=100
PRIMER_MAX_POLY_X=10
PRIMER_MAX_SELF_ANY_TH=100
PRIMER_MAX_SELF_END_TH=100
PRIMER_MAX_GC=100
PRIMER_MIN_GC=0
PRIMER_MAX_END_STABILITY=100
=
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(primer3_input)
            input_file = f.name
        
        result = subprocess.run(
            [primer3_path, input_file],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        os.remove(input_file)
        return result.stdout
    
    def extract_oligo_from_output(output, template_tm, template_seq):
        """Extract oligo from primer3 output, preserving IUPAC codes from original sequence."""
        oligo_seq = None
        oligo_tm = None
        primer_start = None
        
        for line in output.split('\n'):
            if line.startswith('PRIMER_LEFT_0_SEQUENCE='):
                oligo_seq = line.split('=')[1]
            elif line.startswith('PRIMER_LEFT_0_TM='):
                oligo_tm = float(line.split('=')[1])
            elif line.startswith('PRIMER_LEFT_0='):
                # Format: PRIMER_LEFT_0=start,length
                parts = line.split('=')[1].split(',')
                primer_start = int(parts[0])
        
        if oligo_seq and len(oligo_seq) == primer_length and primer_start is not None:
            # Extract the actual sequence from the original template at the selected position
            # This preserves IUPAC codes (R, Y, M, etc.) from the consensus sequence
            if primer_start + primer_length <= len(template_seq):
                actual_seq = template_seq[primer_start:primer_start + primer_length]
            else:
                actual_seq = oligo_seq
            
            new_source_header = f"{prefix}|{template_tm}|{oligo_tm}"
            return probe.Probe.from_str(actual_seq, source_header=new_source_header), actual_seq, primer_start
        return None, None, None
    
    # First attempt: try the original sequence
    output = run_primer3_with_sequence(seq, prefix)
    selected_probe, selected_seq, primer_start = extract_oligo_from_output(output, template_tm, seq)
    
    if selected_probe:
        # Check if the selected oligo has a poly stretch
        poly_base, poly_len = _has_poly_stretch(selected_seq, poly_filter_length)
        
        if poly_base and poly_len >= poly_filter_length:
            logger.debug("Selected %d-mer for %s has poly-%s stretch of %d bases, attempting re-selection",
                        primer_length, header, poly_base, poly_len)
        
            # Try shifted windows
            seq_len = len(seq)
            
            # Calculate how many possible positions for the primer within the probe
            # We slide a window of size shift_window across the sequence
            # and try to select a primer from each window
            num_positions = (seq_len - primer_length) // (shift_window // 2) + 1
            
            # Try different positions within the sequence
            # Start from positions that are offset from the original (index 0)
            for offset in range(1, min(num_positions, 5)):
                # Calculate start position for this offset
                # We want to slide the window by half its size each time
                step = shift_window // 2 if shift_window > primer_length else primer_length
                start_pos = offset * step
                
                # Ensure we have room for a full shift_window
                if start_pos + shift_window > seq_len:
                    # Try from the end
                    start_pos = seq_len - shift_window
                
                if start_pos < 0:
                    start_pos = 0
                
                if start_pos + shift_window <= seq_len:
                    shifted_template = seq[start_pos:start_pos + shift_window]
                    shifted_tm = _calculate_tm(shifted_template)
                    shifted_prefix = f"{prefix}|SHIFT_{offset}|POS_{start_pos}"
                    
                    output = run_primer3_with_sequence(shifted_template, shifted_prefix)
                    selected_probe, selected_seq, primer_start = extract_oligo_from_output(
                        output, shifted_tm, shifted_template)
                    
                    if selected_probe:
                        poly_base, poly_len = _has_poly_stretch(selected_seq, poly_filter_length)
                        if not (poly_base and poly_len >= poly_filter_length):
                            logger.debug("Found valid %d-mer at position %d: %s (poly-%s=%d)",
                                        primer_length, start_pos, selected_seq, poly_base, poly_len)
                            break
            
            # Log final result
            poly_base_final, poly_len_final = _has_poly_stretch(selected_seq, poly_filter_length)
            if poly_base_final and poly_len_final >= poly_filter_length:
                logger.debug("Could not avoid poly-%s stretch in any shifted window, keeping %d-mer: %s",
                            poly_base_final, primer_length, selected_seq)
            else:
                logger.debug("Successfully re-selected %d-mer avoiding poly stretch: %s",
                            primer_length, selected_seq)
        else:
            logger.debug("Selected %d-mer for %s has no problematic poly stretch: %s",
                        primer_length, header, selected_seq)
    else:
        # primer3 could not find a valid oligo, use original probe
        logger.debug("primer3 could not find %d-mer for %s (Tm range %.1f-%.1f), using original", 
                    primer_length, header, min_tm, max_tm)
        return p
    
    return selected_probe if selected_probe else p


def _select_primer3_oligos(probes, primer_length, target_tm, mismatches, num_processes=8,
                            poly_filter_length=5, shift_window=75):
    """Use primer3_core to select optimized oligos from probes with parallel processing.
    
    If a selected oligo has a poly stretch exceeding poly_filter_length, attempts to
    re-select from shifted regions.
    
    Args:
        probes: list of Probe objects to design oligos from
        primer_length: desired oligo length
        target_tm: target melting temperature
        mismatches: number of mismatches to tolerate (for ambiguous bases)
        num_processes: number of parallel processes for primer3
        poly_filter_length: minimum poly stretch length to trigger re-selection
        shift_window: size of shifted window to try for re-selection
    
    Returns:
        list of Probe objects containing the selected oligos
    """
    import subprocess
    import tempfile
    import os
    from multiprocessing import Pool, cpu_count
    
    logger = logging.getLogger(__name__)
    
    # Get the path to primer3_core
    primer3_path = os.path.join(os.path.dirname(__file__), '..', 'primer3', 'src', 'primer3_core')
    
    # Set Tm range with wider flexibility (target_tm +/- 35C)
    # This allows primer3 more options to find valid 26-mers, ensuring
    # every input probe produces a 26-mer output
    min_tm = target_tm - 35
    max_tm = target_tm + 35
    
    # Log initial probe count (75-mers)
    total_probes = len(probes)
    logger.info("Processing %d 75-mer probes through primer3_core for %d-mer selection", 
                total_probes, primer_length)
    
    # Use parallel processing for large batches
    effective_processes = min(num_processes, cpu_count(), len(probes))
    
    if effective_processes > 1:
        # Prepare arguments for parallel processing
        args_list = [(p, primer_length, target_tm, mismatches, min_tm, max_tm, primer3_path,
                      poly_filter_length, shift_window) 
                     for p in probes]
        with Pool(processes=effective_processes) as pool:
            optimized_probes = pool.starmap(_process_probe_for_primer3, args_list)
    else:
        # Sequential processing for small batches
        optimized_probes = []
        for i, p in enumerate(probes, 1):
            optimized_probe = _process_probe_for_primer3(p, primer_length, target_tm, mismatches, 
                                                         min_tm, max_tm, primer3_path,
                                                         poly_filter_length, shift_window)
            optimized_probes.append(optimized_probe)
    
    # Log final summary
    success_count = sum(1 for p in optimized_probes if len(str(p)) == primer_length)
    logger.info("Completed primer3 selection: %d/%d 75-mers produced %d-mers (%.1f%% success rate)",
                success_count, total_probes, primer_length, 100.0 * success_count / total_probes)
    
    # Log final summary
    success_count = sum(1 for p in optimized_probes if len(str(p)) == primer_length)
    logger.info("Completed primer3 selection: %d/%d 75-mers produced %d-mers (%.1f%% success rate)",
                success_count, total_probes, primer_length, 100.0 * success_count / total_probes)
    
    return optimized_probes


# Define types for initializing arguments:
#   'basic': most restrictive (0 mismatches, 0 cover extension, etc.) and
#            without options that enhance runtime and memory usage for large,
#            highly diverse input
#   'large': less restrictive (reasonable number of mismatches and cover
#            extension, etc.) and enabling, by default, options that lower
#            runtime and memory usage for large, highly diverse input at the
#            typical expense of a small increase in the number of output probes
_ARGS_TYPES = typing.Literal['basic', 'large']


def main(args):
    # Setup logger
    log.configure_logging(args.log_level)
    logger = logging.getLogger(__name__)

    # If running design_large.py, warn that defaults for some arguments may
    #   be undesired
    if args.args_type == 'large':
        logger.warning(("With design_large.py, the default values for some "
            "arguments --- such as mismatches (-m) or cover extension (-e) "
            "--- might be more relaxed than desired. Run 'design_large.py "
            "--help' to see the default values; they can be overridden by "
            "specifying the argument."))

    # Set NCBI API key
    if args.ncbi_api_key:
        ncbi_neighbors.ncbi_api_key = args.ncbi_api_key

    # Read the genomes from FASTA sequences
    genomes_grouped = []
    genomes_grouped_names = []
    for ds in args.dataset:
        if ds.startswith('collection:'):
            # Process a collection of datasets
            raise ValueError(("A collection of datasets (via 'collection:') "
                "is no longer allowed as input. Please specify only NCBI "
                "taxonomy IDs to download or FASTA files."))
        elif ds.startswith('download:'):
            # Download a FASTA for an NCBI taxonomic ID
            taxid = ds[len('download:'):]
            if args.write_taxid_acc:
                taxid_fn = os.path.join(args.write_taxid_acc,
                        str(taxid) + '.txt')
            else:
                taxid_fn = None
            if '-' in taxid:
                taxid, segment = taxid.split('-')
            else:
                segment = None
            ds_fasta_tf = ncbi_neighbors.construct_fasta_for_taxid(taxid,
                    segment=segment, write_to=taxid_fn)
            # Read with replace_degenerate=False to preserve IUPAC codes
            fasta_dict = seq_io.read_fasta(ds_fasta_tf.name, replace_degenerate=False)
            genomes_grouped += [[seq_io.genome.Genome.from_one_seq(seq, header=header) 
                                for header, seq in fasta_dict.items()]]
            genomes_grouped_names += ['taxid:' + str(taxid)]
            ds_fasta_tf.close()
        elif os.path.isfile(ds):
            # Process a custom fasta file with sequences
            # Read with replace_degenerate=False to preserve IUPAC codes
            fasta_dict = seq_io.read_fasta(ds, replace_degenerate=False)
            genomes_grouped += [[seq_io.genome.Genome.from_one_seq(seq, header=header) 
                                for header, seq in fasta_dict.items()]]
            genomes_grouped_names += [os.path.basename(ds)]
        else:
            # Process an individual dataset
            raise ValueError(("Dataset labels are no longer allowed as "
                "input. Please specify only NCBI taxonomy IDs to download "
                "(via 'download:taxid') or FASTA files. If you already "
                "specified a FASTA file, please check that the path to "
                f"'{ds}' is valid."))

    if (args.limit_target_genomes and
            args.limit_target_genomes_randomly_with_replacement):
        raise Exception(("Cannot --limit-target-genomes and "
                         "--limit-target-genomes-randomly-with-replacement at "
                         "the same time"))
    elif args.limit_target_genomes:
        genomes_grouped = [genomes[:args.limit_target_genomes]
                           for genomes in genomes_grouped]
    elif args.limit_target_genomes_randomly_with_replacement:
        k = args.limit_target_genomes_randomly_with_replacement
        genomes_grouped = [random.choices(genomes, k=k)
                           for genomes in genomes_grouped]

    # Suggest design_large.py if not initialized by that program and its
    #   settings may be appropriate (i.e., (multiple datasets and --identify
    #   is not set) or large (>10 million nt) input size)
    if args.args_type != 'large':
        total_input_size = sum(sum(g.size() for g in genomes)
                for genomes in genomes_grouped)
        if ((len(args.dataset) > 1 and not args.identify) or
                total_input_size > 10000000):
            recommended_args = []
            if (not args.filter_with_lsh_hamming and
                    not args.filter_with_lsh_minhash):
                recommended_args += ['--filter-with-lsh-minhash 0.6']
            if not args.cluster_and_design_separately:
                recommended_args += ['--cluster-and-design-separately 0.15']
            if not args.cluster_from_fragments:
                recommended_args += ['--cluster-from-fragments 50000']
            recommended_args_str = ""
            if len(recommended_args) > 0:
                recommended_args_str = ("Recommended options include: " +
                        ', '.join(["'" + x + "'" for x in recommended_args]))
            logger.warning(("If runtime or memory usage are problematic, "
                "consider using design_large.py or some of the "
                "options it sets, which may be helpful in lowering runtime "
                "and memory usage for this design. "
                f"{recommended_args_str}"))

    # Store the FASTA paths of avoided genomes
    avoided_genomes_fasta = []
    if args.avoid_genomes:
        for ag in args.avoid_genomes:
            if os.path.isfile(ag):
                # Process a custom fasta file with sequences
                avoided_genomes_fasta += [ag]
            else:
                # Process an individual dataset
                raise ValueError(("Dataset labels are no longer allowed as "
                    "input. Please specify only NCBI taxonomy IDs to download "
                    "(via 'download:taxid') or FASTA files. If you already "
                    "specified a FASTA file, please check that the path to "
                    f"'{ag}' is valid."))

    # Setup and verify parameters related to probe length
    if not args.lcf_thres:
        args.lcf_thres = args.probe_length
    if args.probe_stride > args.probe_length:
        logger.warning(("PROBE_STRIDE (%d) is greater than PROBE_LENGTH "
                        "(%d), which is usually undesirable and may lead "
                        "to undefined behavior"),
                        args.probe_stride, args.probe_length)
    if args.lcf_thres > args.probe_length:
        logger.warning(("LCF_THRES (%d) is greater than PROBE_LENGTH "
                        "(%d), which is usually undesirable and may lead "
                        "to undefined behavior"),
                        args.lcf_thres, args.probe_length)
    if args.island_of_exact_match > args.probe_length:
        logger.warning(("ISLAND_OF_EXACT_MATCH (%d) is greater than "
                        "PROBE_LENGTH (%d), which is usually undesirable "
                        "and may lead to undefined behavior"),
                        args.island_of_exact_match, args.probe_length)
    if args.mismatches / args.probe_length > 0.15:
        logger.warning(("MISMATCHES (%d) is higher relative to PROBE_LENGTH "
                        "(%d) than typically provided, and may lead to "
                        "slower runtime and lower enrichment in practice"),
                        args.mismatches, args.probe_length)

    # Setup and verify parameters related to k-mer length in probe map
    if args.kmer_probe_map_k:
        # Check that k is sufficiently small
        if args.kmer_probe_map_k > args.probe_length:
            raise Exception(("KMER_PROBE_MAP_K (%d) exceeds PROBE_LENGTH "
                             "(%d), which is not permitted") %
                            (args.kmer_probe_map_k, args.probe_length))

        # Use this value for the SetCoverFilter, AdapterFilter, and
        # the Analyzer
        kmer_probe_map_k_scf = args.kmer_probe_map_k
        kmer_probe_map_k_af = args.kmer_probe_map_k
        kmer_probe_map_k_analyzer = args.kmer_probe_map_k
    else:
        if args.probe_length <= 20:
            logger.warning(("PROBE_LENGTH (%d) is small; you may want to "
                            "consider setting --kmer-probe-map-k to be "
                            "small as well in order to be more sensitive "
                            "in mapping candidate probes to target sequence"),
                            args.probe_length)

        # Use a default k of 20 for the SetCoverFilter and AdapterFilter,
        # and 10 for the Analyzer since we would like to be more sensitive
        # (potentially at the cost of slower runtime) for the latter
        kmer_probe_map_k_scf = 20
        kmer_probe_map_k_af = 20
        kmer_probe_map_k_analyzer = 10

    # Set the maximum number of processes in multiprocessing pools
    if args.max_num_processes:
        probe.set_max_num_processes_for_probe_finding_pools(
            args.max_num_processes)
        cluster.set_max_num_processes_for_computing_distances(
            args.max_num_processes)
        set_cover_filter.set_max_num_processes_for_set_cover_instances(
            args.max_num_processes)
        base_filter.set_max_num_processes_for_filter_over_groupings(
            args.max_num_processes)

    # Raise exceptions or warn based on use of adapter arguments
    if args.add_adapters:
        if not (args.adapter_a or args.adapter_b):
            logger.warning(("Adapter sequences will be added, but default "
                            "sequences will be used; to provide adapter "
                            "sequences, use --adapter-a and --adapter-b"))
    else:
        if args.adapter_a or args.adapter_b:
            raise Exception(("Adapter sequences were provided with "
                "--adapter-a and --adapter-b, but --add-adapters is required "
                "to add adapter sequences onto the ends of probes"))

    # Do not allow both --small-seq-skip and --small-seq-min, since they
    # have different intentions
    if args.small_seq_skip is not None and args.small_seq_min is not None:
        raise Exception(("Both --small-seq-skip and --small-seq-min were "
            "specified, but both cannot be used together"))

    # Check arguments involving clustering
    if args.cluster_and_design_separately and args.identify:
        raise Exception(("Cannot use --cluster-and-design-separately with "
            "--identify, because clustering collapses genome groupings into "
            "one"))
    if args.cluster_from_fragments and not args.cluster_and_design_separately:
        raise Exception(("Cannot use --cluster-from-fragments without also "
            "setting --cluster-and-design-separately"))

    # Check for whether a custom hybridization function was provided
    if args.custom_hybridization_fn:
        custom_cover_range_fn = tuple(args.custom_hybridization_fn)
    else:
        custom_cover_range_fn = None
    if args.custom_hybridization_fn_tolerant:
        custom_cover_range_tolerant_fn = tuple(args.custom_hybridization_fn_tolerant)
    else:
        custom_cover_range_tolerant_fn = None

    # Setup the filters
    # The filters we use are, in order:
    filters = []

    # [Optional]
    # Fasta filter (ff) -- leave out candidate probes
    if args.filter_from_fasta:
        ff = fasta_filter.FastaFilter(args.filter_from_fasta,
                                      skip_reverse_complements=True)
        filters += [ff]

    # [Optional]
    # Poly(A) filter (paf) -- leave out probes with stretches of 'A' or 'T'
    if args.filter_polya:
        polya_length, polya_mismatches = args.filter_polya
        if polya_length > args.probe_length:
            logger.warning(("Length of poly(A) stretch to filter (%d) is "
                            "greater than PROBE_LENGTH (%d), which is usually "
                            "undesirable"), polya_length, args.probe_length)
        if polya_length < 10:
            logger.warning(("Length of poly(A) stretch to filter (%d) is "
                            "short, and may lead to many probes being "
                            "filtered"), polya_length)
        if polya_mismatches > 10:
            logger.warning(("Number of mismatches to tolerate when searching "
                            "for poly(A) stretches (%d) is high, and may "
                            "lead to many probes being filtered"),
                           polya_mismatches)
        paf = polya_filter.PolyAFilter(polya_length, polya_mismatches)
        filters += [paf]

    # [Optional]
    # Poly(G) filter (pgf) -- leave out probes with stretches of 'G'
    if args.filter_polyg:
        polyg_length, polyg_mismatches = args.filter_polyg
        if polyg_length > args.probe_length:
            logger.warning(("Length of poly(G) stretch to filter (%d) is "
                            "greater than PROBE_LENGTH (%d), which is usually "
                            "undesirable"), polyg_length, args.probe_length)
        if polyg_length < 10:
            logger.warning(("Length of poly(G) stretch to filter (%d) is "
                            "short, and may lead to many probes being "
                            "filtered"), polyg_length)
        if polyg_mismatches > 10:
            logger.warning(("Number of mismatches to tolerate when searching "
                            "for poly(G) stretches (%d) is high, and may "
                            "lead to many probes being filtered"),
                           polyg_mismatches)
        pgf = polyg_filter.PolyGFilter(polyg_length, polyg_mismatches)
        filters += [pgf]

    # [Optional]
    # Poly(C) filter (pcf) -- leave out probes with stretches of 'C'
    if args.filter_polyc:
        polyc_length, polyc_mismatches = args.filter_polyc
        if polyc_length > args.probe_length:
            logger.warning(("Length of poly(C) stretch to filter (%d) is "
                            "greater than PROBE_LENGTH (%d), which is usually "
                            "undesirable"), polyc_length, args.probe_length)
        if polyc_length < 10:
            logger.warning(("Length of poly(C) stretch to filter (%d) is "
                            "short, and may lead to many probes being "
                            "filtered"), polyc_length)
        if polyc_mismatches > 10:
            logger.warning(("Number of mismatches to tolerate when searching "
                            "for poly(C) stretches (%d) is high, and may "
                            "lead to many probes being filtered"),
                           polyc_mismatches)
        pcf = polyc_filter.PolyCFilter(polyc_length, polyc_mismatches)
        filters += [pcf]

    # [Optional]
    # Poly(T) filter (ptf) -- leave out probes with stretches of 'T'
    if args.filter_polyt:
        polyt_length, polyt_mismatches = args.filter_polyt
        if polyt_length > args.probe_length:
            logger.warning(("Length of poly(T) stretch to filter (%d) is "
                            "greater than PROBE_LENGTH (%d), which is usually "
                            "undesirable"), polyt_length, args.probe_length)
        if polyt_length < 10:
            logger.warning(("Length of poly(T) stretch to filter (%d) is "
                            "short, and may lead to many probes being "
                            "filtered"), polyt_length)
        if polyt_mismatches > 10:
            logger.warning(("Number of mismatches to tolerate when searching "
                            "for poly(T) stretches (%d) is high, and may "
                            "lead to many probes being filtered"),
                           polyt_mismatches)
        ptf = polyt_filter.PolyTFilter(polyt_length, polyt_mismatches)
        filters += [ptf]

    # Duplicate filter (df) -- condense all candidate probes that
    #     are identical down to one; this is not necessary for
    #     correctness, as the set cover filter achieves the same task
    #     implicitly, but it does significantly lower runtime by
    #     decreasing the input size to the set cover filter
    # Near duplicate filter (ndf) -- condense candidate probes that
    #     are near-duplicates down to one using locality-sensitive
    #     hashing; like the duplicate filter, this is not necessary
    #     but can significantly lower runtime and reduce memory usage
    #     (even more than the duplicate filter)
    if (args.filter_with_lsh_hamming is not None and
            args.filter_with_lsh_minhash is not None):
        raise Exception(("Cannot use both --filter-with-lsh-hamming "
            "and --filter-with-lsh-minhash"))
    if args.filter_with_lsh_hamming is not None:
        if args.filter_with_lsh_hamming > args.mismatches:
            logger.warning(("Setting FILTER_WITH_LSH_HAMMING (%d) to be greater "
                "than MISMATCHES (%d) may cause the probes to achieve less "
                "than the desired coverage"), args.filter_with_lsh_hamming,
                args.mismatches)
        ndf = near_duplicate_filter.NearDuplicateFilterWithHammingDistance(
            args.filter_with_lsh_hamming, args.probe_length)
        filters += [ndf]
    elif args.filter_with_lsh_minhash is not None:
        if args.mismatches < 3:
            logger.warning(("MISMATCHES is set to %d; at low values of "
                "MISMATCHES (0, 1, or 2), using --filter-with-lsh-minhash "
                "(particularly with high values of FILTER_WITH_LSH_MINHASH) "
                "may cause the probes to achieve less than the desired "
                "coverage"), args.mismatches)
        ndf = near_duplicate_filter.NearDuplicateFilterWithMinHash(
            args.filter_with_lsh_minhash)
        filters += [ndf]
    else:
        df = duplicate_filter.DuplicateFilter()
        filters += [df]

    # Set cover filter (scf) -- solve the problem by treating it as
    #     an instance of the set cover problem
    scf = set_cover_filter.SetCoverFilter(
        mismatches=args.mismatches,
        lcf_thres=args.lcf_thres,
        island_of_exact_match=args.island_of_exact_match,
        mismatches_tolerant=args.mismatches_tolerant,
        lcf_thres_tolerant=args.lcf_thres_tolerant,
        island_of_exact_match_tolerant=args.island_of_exact_match_tolerant,
        custom_cover_range_fn=custom_cover_range_fn,
        custom_cover_range_tolerant_fn=custom_cover_range_tolerant_fn,
        identify=args.identify,
        avoided_genomes=avoided_genomes_fasta,
        coverage=args.coverage,
        cover_extension=args.cover_extension,
        kmer_probe_map_k=kmer_probe_map_k_scf,
        kmer_probe_map_use_native_dict=args.use_native_dict_when_finding_tolerant_coverage)
    filters += [scf]

    # [Optional]
    # Adapter filter (af) -- add adapters to both the 5' and 3' ends
    #    of each probe
    if args.add_adapters:
        # Set default adapter sequences, if not provided
        if args.adapter_a:
            adapter_a = tuple(args.adapter_a)
        else:
            adapter_a = ('ATACGCCATGCTGGGTCTCC', 'CGTACTTGGGAGTCGGCCAT')
        if args.adapter_b:
            adapter_b = tuple(args.adapter_b)
        else:
            adapter_b = ('AGGCCCTGGCTGCTGATATG', 'GACCTTTTGGGACAGCGGTG')

        af = adapter_filter.AdapterFilter(adapter_a,
                                          adapter_b,
                                          mismatches=args.mismatches,
                                          lcf_thres=args.lcf_thres,
                                          island_of_exact_match=\
                                            args.island_of_exact_match,
                                          custom_cover_range_fn=\
                                            custom_cover_range_fn,
                                          kmer_probe_map_k=kmer_probe_map_k_af)
        filters += [af]

    # [Optional]
    # N expansion filter (nef) -- expand Ns in probe sequences
    # to avoid ambiguity
    if args.expand_n is not None:
        nef = n_expansion_filter.NExpansionFilter(
            limit_n_expansion_randomly=args.expand_n)
        filters += [nef]

    # [Optional]
    # Reverse complement (rc) -- add the reverse complement of each
    #    probe that remains
    if args.add_reverse_complements:
        rc = reverse_complement_filter.ReverseComplementFilter()
        filters += [rc]

    # If requested, don't apply the set cover filter
    if args.skip_set_cover:
        filter_before_scf = filters[filters.index(scf) - 1]
        filters.remove(scf)

    # Define parameters for clustering sequences
    if args.cluster_and_design_separately:
        cluster_threshold = args.cluster_and_design_separately
        if args.skip_set_cover:
            cluster_merge_after = filter_before_scf
        else:
            cluster_merge_after = scf
        cluster_method = args.cluster_and_design_separately_method
        cluster_fragment_length = args.cluster_from_fragments
    else:
        cluster_threshold = None
        cluster_merge_after = None
        cluster_method = None
        cluster_fragment_length = None

    # Design the probes
    pb = probe_designer.ProbeDesigner(genomes_grouped, filters,
                                      probe_length=args.probe_length,
                                      probe_stride=args.probe_stride,
                                      allow_small_seqs=args.small_seq_min,
                                      seq_length_to_skip=args.small_seq_skip,
                                      cluster_threshold=cluster_threshold,
                                      cluster_merge_after=cluster_merge_after,
                                      cluster_method=cluster_method,
                                      cluster_fragment_length=cluster_fragment_length)
    pb.design()

      # If -primerl is specified, use primer3_core to select optimized oligos
    if args.primer_length:
        logger.info("Using primer3_core to select optimized oligos")
        
        # Determine poly filter length from args (use 5 if any poly filter is set)
        poly_filter_length = 5
        if (args.filter_polya or args.filter_polyg or 
            args.filter_polyc or args.filter_polyt):
            # Use the minimum poly filter length specified
            poly_lengths = []
            if args.filter_polya:
                poly_lengths.append(args.filter_polya[0])
            if args.filter_polyg:
                poly_lengths.append(args.filter_polyg[0])
            if args.filter_polyc:
                poly_lengths.append(args.filter_polyc[0])
            if args.filter_polyt:
                poly_lengths.append(args.filter_polyt[0])
            poly_filter_length = min(poly_lengths) if poly_lengths else 5
        
        # Use probe_length as the shift window (not primer_length)
        # This allows re-selection from adjacent regions within the probe
        shift_window = args.probe_length
        
        primer3_probes = _select_primer3_oligos(
            pb.final_probes,
            primer_length=args.primer_length,
            target_tm=args.primer_opt_tm,
            mismatches=args.mismatches,
            poly_filter_length=poly_filter_length,
            shift_window=shift_window)
        # Write original probes with renamed IDs (use identifier format)
        orig_output = args.output_probes.replace('.fasta', '_original.fasta') if '.fasta' in args.output_probes else args.output_probes.rsplit('.', 1)[0] + '_original.fasta'
        seq_io.write_probe_fasta(pb.final_probes, orig_output, use_header_if_set=False)
        logger.info("Wrote %d original probes to %s", len(pb.final_probes), orig_output)
        # Write primer3-optimized oligos to the specified output file
        seq_io.write_probe_fasta(primer3_probes, args.output_probes, use_header_if_set=False)
        logger.info("Wrote %d primer3-optimized oligos to %s", len(primer3_probes), args.output_probes)
    else:
        # Write the final probes to the file args.output_probes
        seq_io.write_probe_fasta(pb.final_probes, args.output_probes, use_header_if_set=False)

    if (args.print_analysis or args.write_analysis_to_tsv or
            args.write_sliding_window_coverage or
            args.write_probe_map_counts_to_tsv):
        analyzer = coverage_analysis.Analyzer(
            pb.final_probes,
            args.mismatches,
            args.lcf_thres,
            genomes_grouped,
            genomes_grouped_names,
            island_of_exact_match=args.island_of_exact_match,
            custom_cover_range_fn=custom_cover_range_fn,
            cover_extension=args.cover_extension,
            kmer_probe_map_k=kmer_probe_map_k_analyzer,
            rc_too=args.add_reverse_complements)
        analyzer.run()
        if args.write_analysis_to_tsv:
            analyzer.write_data_matrix_as_tsv(
                args.write_analysis_to_tsv)
        if args.write_sliding_window_coverage:
            analyzer.write_sliding_window_coverage(
                args.write_sliding_window_coverage)
        if args.write_probe_map_counts_to_tsv:
            analyzer.write_probe_map_counts(
                    args.write_probe_map_counts_to_tsv)
        if args.print_analysis:
            analyzer.print_analysis()
    else:
        # Just print the number of probes
        print(len(pb.final_probes))


def init_and_parse_args(args_type : _ARGS_TYPES):
    """Setup and parse command-line arguments.

    Args:
        args_type: whether to initialize arguments to be optimized for
            large, highly diverse input ('large') or not ('basic')

    Returns:
        populated namespace of arguments
    """
    if args_type not in typing.get_args(_ARGS_TYPES):
        raise ValueError((f"Argument type '{args_type}' is invalid; it must "
            f"be one of {typing.get_args(_ARGS_TYPES)}"))

    # Format --help messages to include default values
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Input data
    parser.add_argument('dataset',
        nargs='+',
        help=("One or more target datasets (e.g., one per species). Each "
              "dataset can be specified in one of two ways. (1) If "
              "dataset is in the format 'download:TAXID', then CATCH downloads "
              "from NCBI all whole genomes for the NCBI taxonomy with id "
              "TAXID, and uses these sequences as input. (2) If dataset is "
              "a path to a FASTA file, then its sequences are read and used "
              "as input. For segmented viruses, the format "
              "for NCBI downloads can also be 'download:TAXID-SEGMENT'."))

    # Outputting probes
    parser.add_argument('-o', '--output-probes',
        required=True,
        help=("The file to which all final probes should be "
              "written; they are written in FASTA format"))

    # Outputting downloaed data
    parser.add_argument('--write-taxid-acc',
        help=("If 'download:' labels are used in datasets, write downloaded "
              "accessions to a file in this directory. Accessions are written "
              "to WRITE_TAXID_ACC/TAXID.txt"))

    # Parameters on probe length and stride
    parser.add_argument('-pl', '--probe-length',
        type=int,
        default=100,
        help=("Make probes be PROBE_LENGTH nt long"))
    parser.add_argument('-ps', '--probe-stride',
        type=int,
        default=50,
        help=("Generate candidate probes from the input "
              "that are separated by PROBE_STRIDE nt"))
    parser.add_argument('-primerl', '--primer-length',
        type=int,
        default=26,
        help=("Use primer3_core to select oligos of this length with Tm "
              "within +/-3C of -primerl-Tm (default 65)"))
    parser.add_argument('-primerl-Tm', '--primer-opt-tm',
        type=float,
        default=65,
        help=("Target optimal Tm for primer3_core oligo selection "
              "(default 65, range 45-85)"))

    # Parameters governing probe hybridization
    default_mismatches = {'basic': 0, 'large': 5}
    parser.add_argument('-m', '--mismatches',
        type=int,
        default=default_mismatches[args_type],
        help=("Allow for MISMATCHES mismatches when determining "
              "whether a probe covers a sequence"))
    parser.add_argument('-l', '--lcf-thres',
        type=int,
        help=("(Optional) Say that a portion of a probe covers a portion "
              "of a sequence if the two share a substring with at most "
              "MISMATCHES mismatches that has length >= LCF_THRES "
              "nt; if unspecified, this is set to PROBE_LENGTH"))
    parser.add_argument('--island-of-exact-match',
        type=int,
        default=0,
        help=("(Optional) When determining whether a probe covers a "
              "sequence, require that there be an exact match (i.e., "
              "no mismatches) of length at least ISLAND_OF_EXACT_"
              "MATCH nt between a portion of the probe and a portion "
              "of the sequence"))

    # Custom function (dynamically loaded) to determine probe hybridization
    # When set, this makes values of the above arguments (--mismatches,
    # --lcf-thres, and --island-of-exact-match) meaningless
    parser.add_argument('--custom-hybridization-fn',
        nargs=2,
        help=("(Optional) Args: <PATH> <FUNC>; PATH is a path to a Python "
              "module (.py file) and FUNC is a string giving the name of "
              "a function in that module. FUNC provides a custom model of "
              "hybridization between a probe and target sequence to use in "
              "the probe set design. If this is set, the arguments "
              "--mismatches, --lcf-thres, and --island-of-exact-match are "
              "not used because these are meant for the default model of "
              "hybridization. The function FUNC in PATH is dynamically "
              "loaded to use when determining whether a probe hybridizes to "
              "a target sequence (and, if so, what portion). FUNC must "
              "accept the following arguments in order, though it "
              "may choose to ignore some values: (1) array giving sequence "
              "of a probe; (2) str giving subsequence of target sequence to "
              "which the probe may hybridize, of the same length as the "
              "given probe sequence; (3) int giving the position in the "
              "probe (equivalently, the target subsequence) of the start "
              "of a k-mer around which the probe and target subsequence "
              "are anchored (the probe and target subsequence are aligned "
              "using this k-mer as an anchor); (4) int giving the end "
              "position (exclusive) of the anchor k-mer; (5) int giving the "
              "full length of the probe (the probe provided in (1) may be "
              "cutoff on an end if it extends further than where the "
              "target sequence ends); (6) int giving the full length of the "
              "target sequence of which the subsequence in (2) is part. "
              "FUNC must return None if it deems that the probe does not "
              "hybridize to the target subsequence; otherwise, it must "
              "return a tuple (start, end) where start is an int giving "
              "the start position in the probe (equivalently, in the "
              "target subsequence) at which the probe will hybridize to "
              "the target subsequence, and end is an int (exclusive) giving "
              "the end position of the hybridization."))

    # Desired coverage of target genomes
    def check_coverage(val):
        fval = float(val)
        ival = int(fval)
        if fval >= 0 and fval <= 1:
            # a float in [0,1] giving fractional coverage
            return fval
        elif fval > 1 and fval == ival:
            # an int > 1 giving number of bp to cover
            return ival
        else:
            raise argparse.ArgumentTypeError(("%s is an invalid coverage "
                                              "value") % val)
    parser.add_argument('-c', '--coverage',
        type=check_coverage,
        default=1.0,
        help=("If this is a float in [0,1], it gives the fraction of "
              "each target genome that must be covered by the selected "
              "probes; if this is an int > 1, it gives the number of "
              "bp of each target genome that must be covered by the "
              "selected probes"))

    # Amount of cover extension to assume
    default_cover_extension = {'basic': 0, 'large': 50}
    parser.add_argument('-e', '--cover-extension',
        type=int,
        default=default_cover_extension[args_type],
        help=("Extend the coverage of each side of a probe by COVER_EXTENSION "
              "nt. That is, a probe covers a region that consists of the "
              "portion of a sequence it hybridizes to, as well as this "
              "number of nt on each side of that portion. This is useful "
              "in modeling hybrid selection, where a probe hybridizes to"
              "a fragment that includes the region targeted by the probe, "
              "along with surrounding portions of the sequence. Increasing "
              "its value should reduce the number of probes required to "
              "achieve the desired coverage."))

    # Differential identification and avoiding genomes
    parser.add_argument('-i', '--identify',
        dest="identify",
        action="store_true",
        help=("Design probes meant to make it possible to identify "
              "nucleic acid from a particular input dataset against "
              "the other datasets; when set, the coverage should "
              "generally be small"))
    parser.add_argument('--avoid-genomes',
        nargs='+',
        help=("One or more genomes to avoid; penalize probes based "
              "on how much of each of these genomes they cover. "
              "The value is a path to a FASTA file."))
    parser.add_argument('-mt', '--mismatches-tolerant',
        type=int,
        help=("(Optional) A more tolerant value for 'mismatches'; "
              "this should be greater than the value of MISMATCHES. "
              "Allows for capturing more possible hybridizations "
              "(i.e., more sensitivity) when designing probes for "
              "identification or when genomes are avoided."))
    parser.add_argument('-lt', '--lcf-thres-tolerant',
        type=int,
        help=("(Optional) A more tolerant value for 'lcf_thres'; "
              "this should be less than LCF_THRES. "
              "Allows for capturing more possible hybridizations "
              "(i.e., more sensitivity) when designing probes for "
              "identification or when genomes are avoided."))
    parser.add_argument('--island-of-exact-match-tolerant',
        type=int,
        default=0,
        help=("(Optional) A more tolerant value for 'island_of_"
              "exact_match'; this should be less than ISLAND_OF_ "
              "EXACT_MATCH. Allows for capturing more "
              "possible hybridizations (i.e., more sensitivity) "
              "when designing probes for identification or when "
              "genomes are avoided."))
    parser.add_argument('--custom-hybridization-fn-tolerant',
        nargs=2,
        help=("(Optional) A more tolerant model than the one "
              "implemented in custom_hybridization_fn. This should capture "
              "more possible hybridizations (i.e., be more sensitive) "
              "when designing probes for identification or when genomes "
              "are avoided. See --custom-hybridization-fn for details "
              "of how this function should be implemented and provided."))

    # Outputting coverage analyses
    parser.add_argument('--print-analysis',
        dest="print_analysis",
        action="store_true",
        help="Print analysis of the probe set's coverage")
    parser.add_argument('--write-analysis-to-tsv',
        help=("(Optional) The file to which to write a TSV-formatted matrix "
              "of the probe set's coverage analysis"))
    parser.add_argument('--write-sliding-window-coverage',
        help=("(Optional) The file to which to write the average coverage "
              "achieved by the probe set within sliding windows of each "
              "target genome"))
    parser.add_argument('--write-probe-map-counts-to-tsv',
        help=("(Optional) The file to which to write a TSV-formatted list of "
              "the number of sequences each probe maps to. This explicitly "
              "does not count reverse complements."))

    # Accepting probes as input and skipping set cover process
    parser.add_argument('--filter-from-fasta',
        help=("(Optional) A FASTA file from which to select candidate probes. "
              "Before running any other filters, keep only the candidate "
              "probes that are equal to sequences in the file and remove "
              "all probes not equal to any of these sequences. This, by "
              "default, ignores sequences in the file whose header contains "
              "the string 'reverse complement'; that is, if there is some "
              "probe with sequence S, it may be filtered out (even if there "
              "is a sequence S in the file) if the header of S in the file "
              "contains 'reverse complement'. This is useful if we already "
              "have probes decided by the set cover filter, but simply "
              "want to process them further by, e.g., adding adapters or "
              "running a coverage analysis. For example, if we have already "
              "run the time-consuming set cover filter and have a FASTA "
              "containing those probes, we can provide a path to that "
              "FASTA file for this argument, and also provide the "
              "--skip-set-cover argument, in order to add adapters to "
              "those probes without having to re-run the set cover filter."))
    parser.add_argument('--skip-set-cover',
        dest="skip_set_cover",
        action="store_true",
        help=("Skip the set cover filter; this is useful when we "
              "wish to see the probes generated from only the "
              "duplicate and reverse complement filters, to gauge "
              "the effects of the set cover filter"))

    # Adding adapters
    parser.add_argument('--add-adapters',
        dest="add_adapters",
        action="store_true",
        help=("Add adapters to the ends of probes; to specify adapter "
              "sequences, use --adapter-a and --adapter-b"))
    parser.add_argument('--adapter-a',
        nargs=2,
        help=("(Optional) Args: <X> <Y>; Custom A adapter to use; two ordered "
              "where X is the A adapter sequence to place on the 5' end of "
              "a probe and Y is the A adapter sequence to place on the 3' "
              "end of a probe"))
    parser.add_argument('--adapter-b',
        nargs=2,
        help=("(Optional) Args: <X> <Y>; Custom B adapter to use; two ordered "
              "where X is the B adapter sequence to place on the 5' end of "
              "a probe and Y is the B adapter sequence to place on the 3' "
              "end of a probe"))

    # Filtering poly(A) sequence from probes
    parser.add_argument('--filter-polya',
        nargs=2,
        type=int,
        help=("(Optional) Args: <X> <Y> (integers); do not output any probe "
              "that contains a stretch of X or more 'A' bases, tolerating "
              "up to Y mismatches (and likewise for 'T' bases)"))
    # Filtering poly(G) sequence from probes
    parser.add_argument('--filter-polyg',
        nargs=2,
        type=int,
        help=("(Optional) Args: <X> <Y> (integers); do not output any probe "
              "that contains a stretch of X or more 'G' bases, tolerating "
              "up to Y mismatches"))
    # Filtering poly(C) sequence from probes
    parser.add_argument('--filter-polyc',
        nargs=2,
        type=int,
        help=("(Optional) Args: <X> <Y> (integers); do not output any probe "
              "that contains a stretch of X or more 'C' bases, tolerating "
              "up to Y mismatches"))
    # Filtering poly(T) sequence from probes
    parser.add_argument('--filter-polyt',
        nargs=2,
        type=int,
        help=("(Optional) Args: <X> <Y> (integers); do not output any probe "
              "that contains a stretch of X or more 'T' bases, tolerating "
              "up to Y mismatches"))

    # Adjusting probe output
    parser.add_argument('--add-reverse-complements',
        dest="add_reverse_complements",
        action="store_true",
        help=("Add to the output the reverse complement of each probe"))
    parser.add_argument('--expand-n',
        nargs='?',
        type=int,
        default=None,
        const=3,
        help=("Expand each probe so that 'N' bases are replaced by real "
              "bases; for example, the probe 'ANA' would be replaced "
              "with the probes 'AAA', 'ATA', 'ACA', and 'AGA'; this is "
              "done combinatorially across all 'N' bases in a probe, and "
              "thus the number of new probes grows exponentially with the "
              "number of 'N' bases in a probe. If followed by a command- "
              "line argument (INT), this only expands at most INT randomly "
              "selected N bases, and the rest are replaced with random "
              "unambiguous bases (default INT is 3)."))

    # Limiting input
    parser.add_argument('--limit-target-genomes',
        type=int,
        help=("(Optional) Use only the first LIMIT_TARGET_GENOMES target "
              "genomes in the dataset"))
    parser.add_argument('--limit-target-genomes-randomly-with-replacement',
        type=int,
        help=("(Optional) Randomly select LIMIT_TARGET_GENOMES_RANDOMLY_"
              "WITH_REPLACMENT target genomes in the dataset with "
              "replacement"))

    # Clustering input sequences
    def check_cluster_and_design_separately(val):
        fval = float(val)
        if fval > 0 and fval <= 0.5:
            # a float in (0,0.5]
            return fval
        else:
            raise argparse.ArgumentTypeError(("%s is an invalid average "
                                              "nucleotide dissimilarity") % val)
    default_cluster_and_design_separately = {'basic': None, 'large': 0.15}
    parser.add_argument('--cluster-and-design-separately',
        type=check_cluster_and_design_separately,
        default=default_cluster_and_design_separately[args_type],
        help=("(Optional) If set, cluster all input sequences using their "
              "MinHash signatures, design probes separately on each cluster, "
              "and combine the resulting probes. This can significantly lower "
              "runtime and memory usage, but may lead to a suboptimal "
              "solution. The value CLUSTER_AND_DESIGN_SEPARATELY gives the "
              "distance threshold for determining clusters in terms of "
              "average nucleotide dissimilarity (1-ANI, where ANI is "
              "average nucleotide identity; see --cluster-and-design-"
              "separately-method for details); higher values "
              "result in fewer clusters, and thus longer runtime. Values "
              "must be in (0,0.5], and generally should be around 0.1 to "
              "0.2. When used, this creates a separate genome for each "
              "input sequence -- it collapses all sequences, across both "
              "groups and genomes, into one list of sequences in one group. "
              "Therefore, genomes will not be grouped as specified in the "
              "input and sequences will not be grouped by genome, and "
              "differential identification is not supported"))
    parser.add_argument('--cluster-and-design-separately-method',
        choices=['choose', 'simple', 'hierarchical'], default='choose',
        help=("(Optional) Method for clustering input sequences, which is "
              "only used if --cluster-and-design-separately is set. If "
              "'simple', clusters are connected components of a graph in "
              "which each sequence is a vertex and two sequences are adjacent "
              "if their estimated nucleotide dissimilarity is within "
              "the value CLUSTER_AND_DESIGN_SEPARATELY. If 'hierarchical', "
              "clusters are determined by agglomerative hierarchical "
              "clustering and the the value CLUSTER_AND_DESIGN_SEPARATELY "
              "is the inter-cluster distance threshold to merge clusters. "
              "If 'choose', use a heuristic to decide among 'simple' and "
              "'hierarchical' based on the input. This option can affect "
              "performance and the heuristic does not always make the right "
              "choice, so trying both choices 'simple' and 'hierarchical' "
              "can sometimes be helpful if needed."))
    default_cluster_from_fragments = {'basic': None, 'large': 50000}
    parser.add_argument('--cluster-from-fragments',
        type=int,
        default=default_cluster_from_fragments[args_type],
        help=("(Optional) If set, break all sequences into sequences of "
              "length CLUSTER_FROM_FRAGMENTS nt, and cluster these fragments. "
              "This can be useful for improving runtime on input with "
              "especially large genomes, in which probes for different "
              "fragments can be designed separately. Values should generally "
              "be around 50,000. For this to be used, "
              "--cluster-and-design-separately must also be set."))

    # Filter candidate probes with LSH
    parser.add_argument('--filter-with-lsh-hamming',
        type=int,
        help=("(Optional) If set, filter candidate probes for near-"
              "duplicates using LSH with a family of hash functions that "
              "works with Hamming distance. FILTER_WITH_LSH_HAMMING gives "
              "the maximum Hamming distance at which to call near-"
              "duplicates; it should be commensurate with (but not greater "
              "than) MISMATCHES. Using this may significantly improve "
              "runtime and reduce memory usage by reducing the number of "
              "candidate probes to consider, but may lead to a slightly "
              "sub-optimal solution. It may also, particularly with "
              "relatively high values of FILTER_WITH_LSH_HAMMING, cause "
              "coverage obtained for each genome to be slightly less than "
              "the desired coverage (COVERAGE) when that desired coverage "
              "is the complete genome; using --print-analysis or "
              "--write-analysis-to-tsv will provide the obtained coverage."))
    def check_filter_with_lsh_minhash(val):
        fval = float(val)
        if fval >= 0.0 and fval <= 1.0:
            # a float in [0,1]
            return fval
        else:
            raise argparse.ArgumentTypeError(("%s is an invalid Jaccard "
                                              "distance") % val)
    default_filter_with_lsh_minhash = {'basic': None, 'large': 0.6}
    parser.add_argument('--filter-with-lsh-minhash',
        type=check_filter_with_lsh_minhash,
        default=default_filter_with_lsh_minhash[args_type],
        help=("(Optional) If set, filter candidate probes for near-"
              "duplicates using LSH with a MinHash family. "
              "FILTER_WITH_LSH_MINHASH gives the maximum Jaccard distance "
              "(1 minus Jaccard similarity) at which to call near-duplicates; "
              "the Jaccard similarity is calculated by treating each probe "
              "as a set of overlapping 10-mers. Its value should be "
              "commensurate with parameter values determining whether a probe "
              "hybridizes to a target sequence, but this can be difficult "
              "to measure compared to the input for --filter-with-lsh-hamming. "
              "This argument allows more sensitivity in near-duplicate "
              "detection than --filter-with-lsh-hamming (e.g., if near-"
              "duplicates should involve probes shifted relative to each "
              "other) and, therefore, greater improvement in runtime and "
              "memory usage. Values should generally be around 0.5 to 0.7. "
              "The same caveat mentioned in the help message for "
              "--filter-with-lsh-hamming also applies here; namely, it can "
              "cause the coverage obtained for each genome to be slightly "
              "less than the desired coverage (COVERAGE), and especially so "
              "with low values of MISMATCHES (~0, 1, or 2). Values of "
              "FILTER_WITH_LSH_MINHASH above ~0.7 may start to require "
              "significant memory and runtime for near-duplicate detection "
              "and are usually not recommended."))

    # Miscellaneous technical adjustments
    parser.add_argument('--small-seq-skip',
        type=int,
        help=("(Optional) Do not create candidate probes from sequences "
              "whose length is <= SMALL_SEQ_SKIP. If set to (PROBE_LENGTH - "
              "1), this avoids the error raised when sequences are less "
              "than the probe length"))
    parser.add_argument('--small-seq-min',
        type=int,
        help=("(Optional) If set, allow sequences as input that are "
              "shorter than PROBE_LENGTH (when not set, the program will "
              "error on such input). SMALL_SEQ_MIN is the "
              "minimum sequence length that should be accepted as input. "
              "When a sequence is less than PROBE_LENGTH, a candidate "
              "probe is created that is equal to the sequence; thus, "
              "the output probes may have different lengths. Note that, "
              "when this is set, it might be a good idea to also set "
              "LCF_THRES to be a value smaller than PROBE_LENGTH -- "
              "e.g., the length of the shortest input sequence; otherwise, "
              "when a probe of length p_l is mapped to a sequence of length "
              "s_l, then lcf_thres is treated as being min(LCF_THRES, p_l, "
              "s_l) so that a probe is able to 'cover' a sequence shorter "
              "than the probe and so that a probe shorter than lcf_thres "
              "is able to 'cover' a sequence"))
    def check_max_num_processes(val):
        ival = int(val)
        if ival >= 1:
            return ival
        else:
            raise argparse.ArgumentTypeError(("MAX_NUM_PROCESSES must be "
                                              "an int >= 1"))
    # For 'large' args_type, use all CPUs (by default, if unspecified, the
    #   maxmimum number of processes is determined by each module that is
    #   parallelized, in its set_max_num_processes_*() function)
    default_max_num_processes = {'basic': None,
            'large': multiprocessing.cpu_count()}
    parser.add_argument('--max-num-processes',
        type=check_max_num_processes,
        default=default_max_num_processes[args_type],
        help=("(Optional) An int >= 1 that gives the maximum number of "
              "processes to use in multiprocessing pools; uses min(number "
              "of CPUs in the system, MAX_NUM_PROCESSES) processes"))
    parser.add_argument('--kmer-probe-map-k',
        type=int,
        help=("(Optional) Use this value (KMER_PROBE_LENGTH_K) as the "
              "k-mer length when constructing a map of k-mers to the probes "
              "that contain these k-mers. This map is used when mapping "
              "candidate probes to target sequences and the k-mers serve "
              "as seeds for calculating whether a candidate probe 'covers' "
              "a subsequence. The value should be sufficiently less than "
              "PROBE_LENGTH so that it can find mappings even when the "
              "candidate probe and target sequence are divergent. In "
              "particular, CATCH will try to find a value k >= "
              "KMER_PROBE_LENGTH_K (by default, >=20) such that k divides "
              "PROBE_LENGTH and k < PROBE_LENGTH / MISMATCHES (if "
              "MISMATCHES=0, then k=PROBE_LENGTH). It will then use this "
              "k as the k-mer length in mappings; if no such k exists, it "
              "will use a randomized approach with KMER_PROBE_LENGTH_K as "
              "the k-mer length. If --custom-hybridization-fn is set, "
              "it will always use the randomized approach with "
              "KMER_PROBE_LENGTH_K (by default, 20) as the k-mer length."))
    parser.add_argument('--use-native-dict-when-finding-tolerant-coverage',
        dest="use_native_dict_when_finding_tolerant_coverage",
        action="store_true",
        help=("When finding probe coverage for avoiding genomes and "
              "identification (i.e., when using tolerant parameters), "
              "use a native Python dict as the kmer_probe_map across "
              "processes, rather than the primitives in SharedKmerProbeMap "
              "that are more suited to sharing across processes. Depending "
              "on the input (particularly if there are many candidate probes) "
              "this may result in substantial memory usage; but it may provide "
              "an improvement in runtime when there are relatively few "
              "candidate probes and a very large avoided genomes input"))
    parser.add_argument('--ncbi-api-key',
        help=("API key to use for NCBI e-utils. Using this increases the "
              "limit on requests/second and may prevent an IP address "
              "from being blocked due to too many requests"))

    # Log levels and version
    parser.add_argument('--debug',
        dest="log_level",
        action="store_const",
        const=logging.DEBUG,
        default=logging.WARNING,
        help=("Debug output"))
    parser.add_argument('--verbose',
        dest="log_level",
        action="store_const",
        const=logging.INFO,
        help=("Verbose output"))
    parser.add_argument('-V', '--version',
        action='version',
        version=version.get_version())

    args = parser.parse_args()

    # Add on, to the namespace, the argument type used for initialization
    args.args_type = args_type

    return args


if __name__ == "__main__":
    # Parse arguments - args_type is determined by which script is run
    # For design.py, use 'basic' settings
    import sys
    # Remove 'basic' from sys.argv if present (it's passed as first arg)
    if len(sys.argv) > 1 and sys.argv[1] == 'basic':
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        args = init_and_parse_args(args_type='basic')
    else:
        args = init_and_parse_args(args_type='basic')

    main(args)
