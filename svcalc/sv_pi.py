#!/usr/bin/env python3

"""

This module contains functions for calculating the
sv pi statistic in sliding windows across the genome.

It contains one class, VarCentric, which is a subclass
of VcfSvPositions. This class is only used to reorient
the positions dictionary to be chromosome -> variant_id.

The main function in this module is calculate_sv_pi_full,
which calculates the sv pi statistic for sliding windows,
using the VarCentric class and the functions get_complex_pairwise,
and calculate_sv_pi_windows.

"""
import os
from itertools import combinations
from typing import List, Dict, Tuple
import numpy
import egglib
import bdsg
from .sv_utilities import check_sv, apply_jukes_cantor, \
                          verbose_function, reformat_vcf_egglib, \
                          split_vcf
from .sv_to_bed import VcfSvPositions

class VarCentric(VcfSvPositions):

    """

    This class is a subclass of VcfSvPositions.
    It is simply used to reorient the positions dictionary
    to be chromosome -> variant_id.

    Parameters:

    (Same as VcfSvPositions)

    Returns:

    None

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parse_vcf()
        self.reoriented_positions = dict()

    # Reorients the dictionary of positions to be chromosome -> varid -> individual -> [start, end].
    def reorient_positions_dict(self):
        """

        This takes the positions dictionary built
        in the VcfSvPositions class and reorients
        it to be chromosome -> variant_id -> individual -> [start, end].

        """

        for individual, individual_dict in self.vcf_positions_dict.items():

            for chromosome, positions in individual_dict.items():

                for each in positions:
                    (start, end, node_id) = each

                    if chromosome not in self.reoriented_positions:
                        self.reoriented_positions[chromosome] = dict()

                    if node_id not in self.reoriented_positions[chromosome]:
                        self.reoriented_positions[chromosome][node_id] = dict()

                    self.reoriented_positions[chromosome][node_id][individual] = [start, end]

def get_complex_pairwise(vcfs: list, verbose=True):

    """

    This function gets the positions of SVs for all pairwise comparisons.

    Parameters:

    vcfs:    A list of VCF files to be compared. Must be of the format individual1_individual2.vcf.
    verbose: If True, print messages to stderr to show progress.

    Returns:

    sv_positions_pairwise: A dictionary of the form
    chromosome -> individual_a -> individual_b -> set of positions that are SVs.

    """

    verbose_function(verbose, "Getting positions of SVs for all pairwise comparisons.")
    sv_positions_pairwise = dict()

    for vcf in vcfs:
        verbose_function(verbose, "Opening VCF file " + vcf + " with egglib.")
        (indiv_a, indiv_b) = os.path.basename(vcf).strip().split(".")[0].split("_")[0:2]
        verbose_function(verbose, "Getting SV positions for " + indiv_a + " and " + indiv_b + ".")
        vcf_egg = egglib.io.VCF(vcf)

        while vcf_egg.read():
            chromosome = vcf_egg.get_chrom()
            position = vcf_egg.get_pos()
            alleles = vcf_egg.get_alleles()
            is_an_sv = check_sv(alleles[0], ",".join(alleles[1]))

            if is_an_sv:

                if chromosome not in sv_positions_pairwise:
                    sv_positions_pairwise[chromosome] = dict()

                if indiv_a not in sv_positions_pairwise[chromosome]:
                    sv_positions_pairwise[chromosome][indiv_a] = dict()

                if indiv_b not in sv_positions_pairwise[chromosome][indiv_a]:
                    sv_positions_pairwise[chromosome][indiv_a][indiv_b] = set()

                sv_positions_pairwise[chromosome][indiv_a][indiv_b].add(position)

    return sv_positions_pairwise

def calculate_sv_pi_windows(genotypes_dict: Dict[str, Dict[str, List[str]]],
                            vcf: str,
                            reoriented_positions: Dict[str, Dict[str, Dict[str, Tuple[int, int]]]],
                            sv_positions_pairwise: Dict[str, Dict[str, Dict[str, set]]],
                            complex_variants: set,
                            window_size=100000,
                            step_size=10000,
                            verbose=True):
    """

    This function calculates the sv pi statistic for sliding windows across the genome.

    Parameters:

    genotypes_dict:        Genotypes dictionary of the form
                           chromosome -> variant_id -> genotype_list.
                           Needed because the VCF file is not in the correct
                           format for egglib.

    vcf:                   The name of the VCF file to be reformatted and
                           passed to egglib.

    reoriented_positions:  Returned by the reorient_positions_dict method of
                           the VarCentric class. This is a dictionary of the form
                           chromosome -> variant_id -> individual -> [start, end].

    sv_positions_pairwise: Returned by the get_complex_pairwise function. This
                           is a dictionary of the form
                           chromosome -> individual_a -> individual_b ->
                           set of positions that are SVs.
                           All sites in individual_a that are SVs between individual_a
                           and individual_b are recorded.
                           The idea of this function is to calculate the average number
                           of pairwise SVs for complex variants using this dictionary,
                           and for simple variants using a simple combination of genotypes.

    complex_variants:      A set of variant IDs that are complex variants. These are
                           variants that have more than one allele.

    verbose:               If True, print messages to stderr to show progress.

    Returns:

    None

    """

    verbose_function(verbose, ("Doing sliding window analysis. Window size = "
                               + str(window_size) + "; step size = " + str(step_size) + "."))
    verbose_function(verbose, "Opening VCF file with egglib.")

    vcf_egg = egglib.io.VcfParser(vcf)

    if not os.path.exists(vcf+"i"):
        verbose_function(verbose, "Index for VCF file not found. Creating index.")
        egglib.io.make_vcf_index(vcf)

    verbose_function(verbose, "Loading index" + vcf + "i" + " for VCF file " + vcf + ".")
    vcf_egg.load_index(vcf+"i")

    verbose_function(verbose, "Calculating sv pi.")

    for window in vcf_egg.slider(window_size,
                                 step_size,
                                 start=0,
                                 allow_indel=True,
                                 allow_custom=True):

        # This is the chromosome that the window is on.
        chromosome = window.chromosome

        # The structural variant count starts at zero.
        # This is the average number of pairwise SVs in the window.
        variant_count = 0

        # Window size is the window size minus the number of variant sites in the window.
        window_size_corrected = window_size - window.num_sites

        # This is a dictionary of the form individual -> window_size.
        # We need to calculate the average window size across all individuals.
        all_window_sizes = dict()

        verbose_function(verbose,
                         (str(window.num_sites) + " sites in window " +
                          chromosome + " " + str(window.bounds[0]) +
                          " to " + str(window.bounds[1])))

        for site in window:

            # Start by assuming the site is not an SV.
            svar = False

            # Go to the site in the VCF file.
            vcf_egg.goto(chromosome, int(site.position))

            # Read the site.
            vcf_egg.readline()

            # Get the variant at the site.
            variant = vcf_egg.get_variant()

            # Get the positions of the variant in the reoriented positions dictionary.
            # For each individual, we have a start and end position for the variant.
            this_chromosome_dict = reoriented_positions[chromosome]
            this_variant_positions = this_chromosome_dict[variant.ID[0]]

            # For all individuals' positions.
            for individual, (start, end) in this_variant_positions.items():

                # If the individual has not had its window size recorded,
                # start with the window size minus the number of variant sites.
                if individual not in all_window_sizes:
                    all_window_sizes[individual] = window_size_corrected

                # For this individual, add the length of the variant to the window size.
                all_window_sizes[individual] += (end - start)

                # If the variant length for any individual is greater than or
                # equal to 50, it is an SV.
                if (end - start) >= 50:
                    svar = True

            # If the variant is an SV, we need to calculate the average
            # number of pairwise SVs for the site.
            # Start with a position count, i.e. number of pairwise SVs for the site, of zero.
            position_count = 0

            # Only if the site is an sv.
            if svar:

                # If the variant is complex, we need to calculate the average number
                # of pairwise SVs for the site from the sv_positions_pairwise dictionary.
                if variant.ID[0] in complex_variants:
                    # This is the dictionary of pairwise SVs for the chromosome,
                    # across all genotype comparisons.
                    this_pairs = sv_positions_pairwise[chromosome]
                    genos = ['0'] + genotypes_dict[chromosome][variant.ID[0]]

                    # For each comparison, we need to calculate the number
                    # of SVs in the complex variant.
                    for indiv_a, pairwise_comparisons in this_pairs.items():

                        for indiv_b in pairwise_comparisons.keys():
                            # So, we get the variant positions for the first individual.
                            # Ordering of indiv_a and indiv_b is arbitrary.
                            # All individuals will be in this_variant_positions, so for
                            # whichever was treated as reference in the comparison, we
                            # just count the number of SV sites
                            # (see function get_complex_pairwise).
                            indiv_a_positions = this_pairs[indiv_a][indiv_b]

                            # Get start and end coordinates for the variant at this site
                            # for indiv_a.
                            if indiv_a in this_variant_positions:
                                (start, end) = this_variant_positions[indiv_a]
                                # Increment the position count by the number of SV sites in the
                                # complex variant for this comparison.
                                position_count += len(set(filter(lambda pos: start <= pos <= end,
                                                                 indiv_a_positions)))

                    # Divide the position count by the number of comparisons
                    # to get the average number of pairwise SVs for the site.
                    position_count = position_count / ((len(genos)*(len(genos)-1))/2)

                else:
                    # If it is not a complex variant, we just need to calculate the
                    # average number of pairwise SVs for the site from the genotypes dictionary.
                    genos = ['0'] + genotypes_dict[chromosome][variant.ID[0]]

                    pairs = list(combinations(genos, 2)) # All pairwise comparisons.
                                                         # Not including self-comparisons.

                    position_count = sum(1 for a, b in pairs if a != b) / \
                                     ((len(genos)*(len(genos)-1))/2) # Position count is just the
                                                             # average number of pairwise
                                                             # SVs for the site.

            variant_count += position_count # Increment the total variant count.
                                            # for the window by the position count.

        # Adjust the window size to be the average window size across all individuals.
        if all_window_sizes:
            window_size_corrected = numpy.mean(list(all_window_sizes.values())) / 1000

        else:
            verbose_function(verbose,
                             ("No window sites found for window " +
                              chromosome + " " +
                              str(window.bounds[0]) +
                              " to " + str(window.bounds[1]) + "."))
            window_size_corrected = window_size / 1000

        # Then we divide the variant count by the average window size,
        # in kb, to get the average number of pairwise SVs for the window.
        if variant_count > 0:
            variant_count_kb = variant_count / window_size_corrected

        else:
            variant_count_kb = 0

        # Then print out the calculations across windows.
        output_values = [variant_count_kb, window_size_corrected,
                         chromosome, window.bounds[0], window.bounds[1]]
        output_values_str = list(map(str, output_values))
        print("\t".join(output_values_str))

def calculate_sv_pi_full(vcf_file: str,
                         gbz_file: str,
                         vcf_files: List[str],
                         temp_dir: str,
                         reference_name: str,
                         verbose=False,
                         threads=1,
                         window_size=100000,
                         step_size=10000):

    """

    This function uses the VarCentric class to calculate the sv pi statistic for
    sliding windows across the genome.

    Parameters:

    vcf_file:       The name of the VCF file to be reformatted and passed to egglib.
    gbz_file:       The name of the GBZ file to be used to get coordinates.
    vcf_files:      The VCF files from the all vs all cactus comparison.
    temp_dir:       The temporary directory to store files.
    reference_name: The name of the reference individual that the VCF was deconstructed against.
                    Should match one of the paths in the GBZ.
    verbose:        If True, print messages to stderr to show progress.
    threads:        The number of threads to use for conversion of the gbz with vg.
    window_size:    The size of the sliding window.
    step_size:      The step size of the sliding window.

    Returns:

    None

    """

    # Reformat the VCF file to be compatible with egglib.
    reformat_vcf_egglib(vcf_file, temp_dir)
    vcf_file_for_egg = temp_dir + "/temp_egg.vcf"

    # Get the dictionary of variant positions and reorient them.
    vcf_dict = VarCentric(vcf_file,
                          reference_name,
                          gbz_file,
                          temp_dir,
                          threads=threads,
                          svs=False,
                          verbose=verbose)

    vcf_dict.reorient_positions_dict()
    chromosomes = list(vcf_dict.reoriented_positions.keys())

    sv_positions_pairwise = get_complex_pairwise(vcf_files, verbose=verbose)

    for chromosome in chromosomes:
        split_vcf(vcf_file_for_egg, chromosome, temp_dir)
        vcf_file_for_egg_current = temp_dir+"/"+chromosome+".vcf"
        calculate_sv_pi_windows(vcf_dict.genotypes_dict,
                                vcf_file_for_egg_current,
                                vcf_dict.reoriented_positions,
                                sv_positions_pairwise,
                                vcf_dict.complex_variants,
                                window_size,
                                step_size,
                                verbose)
