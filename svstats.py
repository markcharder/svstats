#!/usr/bin/env python3

"""

This is the main script for the svstats package.
It is a command line tool with multiple functions,
but was designed originally to calculate svpi.

This script imports from the svcalc module, which
is a collection of scripts for working with structural
variants.

Returns:

None

"""

import argparse
import os
import shutil
from svcalc import sv_to_bed, sv_pi, sv_utilities, PACKAGE_VERSION
from svcalc import sv_gene_centric

PARSER = argparse.ArgumentParser(prog="svstats")

PARSER.add_argument("--version",
                    action="version",
                    version=PACKAGE_VERSION)
PARSER.add_argument("--temp-dir",
                    help="Temporary directory. (default = svstats/tmp/)",
                    default=os.getcwd()+"/tmp")
PARSER.add_argument("--threads",
                    help="Threads to use for vg. (default = 1)",
                    default=1)
PARSER.add_argument("--verbose",
                    action="store_true",
                    help="Prints more information during run.",
                    default=False)
PARSER.add_argument("--vg-path",
                    help="Path to vg executable, excluding the executable file itself. \
                          Should be called 'vg'. (default = svstats/bin)",
                    default=os.path.abspath(__file__)+"/bin")

SUBPARSERS = PARSER.add_subparsers(dest="cmd")
SUBPARSERS.required = True

PARSER_A = SUBPARSERS.add_parser("svpi", help="Calculates svpi in sliding windows.")
PARSER_A.add_argument("--step",
                      help="Step size of sliding window in bp. \
                            (default = 5000)",
                      type=int,
                      default=5000)
PARSER_A.add_argument("--window",
                      help="Window size of sliding window in bp. \
                            (default = 50000)",
                      type=int,
                      default=50000)
PARSER_A.add_argument("reference_name",
                      help="The name of the reference sample.")
PARSER_A.add_argument("vcf",
                      help="A VCF file to use if a gbz is not available.")
PARSER_A.add_argument("gbz",
                      help="A GBZ file to get the coordinates of variants in individuals. \
                            See notes on running cactus-pangenome first.")
PARSER_A.add_argument("vcf_files", nargs="+",
                      help="A list of VCF files corresponding to the individuals in the graph.")

PARSER_B = SUBPARSERS.add_parser("svcoords", help="For each structural variant in a VCF, \
                                                   prints a BED file with coordinates \
                                                   for each individual.")
PARSER_B.add_argument("reference_name",
                      help="The name of the reference individual that the \
                            VCF was deconstructed against. \
                            Should match one of the paths in the GBZ.")
PARSER_B.add_argument("vcf",
                      help="A VCF file with node IDs in the third column. \
                            Should correspond to the gbz file provided.")
PARSER_B.add_argument("gbz",
                      help="A GBZ file with all individuals in the VCF as reference sense paths. \
                            Names of individuals should match those in the VCF.")
PARSER_B.add_argument("--sv-size",
                      help="The minimum size of structural variants to include.",
                      default=50)
PARSER_B.add_argument("--consider-others",
                      help="Consider non-structural-variants in the output. \
                            Though originally used only for SVs, you can get the coordinates \
                            for any variant type.",
                      action="store_true",
                      default=False)

PARSER_C = SUBPARSERS.add_parser("randcoords",
                                 help="For an arbitrary BED file, \
                                       get a set of randomised coordinates.")
PARSER_C.add_argument("--n-shuffles",
                      help="The number of randomisations.",
                      default=1000)
PARSER_C.add_argument("bed",
                      help="Any BED file. Only the first three columns are used.")
PARSER_C.add_argument("chromosome_lengths",
                      help="A file with chromosome lengths in it.")

PARSER_D = SUBPARSERS.add_parser("getchromlens",
                                 help="For an arbitrary FASTA file, \
                                       get a TSV with chromosome and length.")
PARSER_D.add_argument("fasta",
                      help="Any FASTA file.")

PARSER_E = SUBPARSERS.add_parser("plink2ldhat",
                                 help="Converts PLINK files to LDhat format.")
PARSER_E.add_argument("--ldhot",
                      action="store_true",
                      default=False,
                      help="Convert to ldhot format instead of ldhat.")
PARSER_E.add_argument("ped",
                      help="A plink PED file.")
PARSER_E.add_argument("map",
                      help="A plink MAP file.")

GROUP_E = PARSER_E.add_mutually_exclusive_group(required=True)
GROUP_E.add_argument("--chromosome-lengths", help="A file with chromosome lengths in it.")
GROUP_E.add_argument("--fasta", help="A FASTA file for creating chromosome lengths.")

PARSER_F = SUBPARSERS.add_parser("watfsites",
                                 help="Calculates a finite sites version of the \
                                       Watterson estimator. Requires ldhat sites \
                                       and locs files.")

PARSER_F.add_argument("sites",
                      help="An LDhat sites file.")
PARSER_F.add_argument("locs",
                      help="An LDhat locs file matching the sites file.")

PARSER_G = SUBPARSERS.add_parser("findcblocks",
                                 help="Finds contiguous blocks of missing or present genes \
                                       in a gene presence / absence matrix produced by \
                                       pangene.")
PARSER_G.add_argument("bubbles",
                      help="The bubbles.txt file output by pangene - follow pangene's GitHub.")

OPTIONS = PARSER.parse_args()

def main():
    """

    This is the main function.
    It tests which subparser \
    is set and runs the appropriate functions.

    Returns:

    None

    """
    if not os.path.exists(OPTIONS.temp_dir):
        os.makedirs(OPTIONS.temp_dir)

    if OPTIONS.cmd == "svpi":
        OPTIONS.gbz = os.path.abspath(OPTIONS.gbz)
        OPTIONS.vcf = os.path.abspath(OPTIONS.vcf)
        sv_pi.calculate_sv_pi_full(OPTIONS.vcf,
                                   OPTIONS.gbz,
                                   OPTIONS.vcf_files,
                                   OPTIONS.temp_dir,
                                   OPTIONS.reference_name,
                                   OPTIONS.verbose,
                                   OPTIONS.threads,
                                   OPTIONS.window,
                                   OPTIONS.step)

    elif OPTIONS.cmd == "svcoords":
        OPTIONS.gbz = os.path.abspath(OPTIONS.gbz)
        OPTIONS.vcf = os.path.abspath(OPTIONS.vcf)
        vcf_sv_positions = sv_to_bed.VcfSvPositions(OPTIONS.vcf, OPTIONS.reference_name,
                                                    OPTIONS.gbz, OPTIONS.temp_dir,
                                                    threads=OPTIONS.threads,
                                                    sv_size=OPTIONS.sv_size,
                                                    consider_others=OPTIONS.consider_others)
        vcf_sv_positions.parse_vcf()
        vcf_sv_positions.print_bed()

    elif OPTIONS.cmd == "randcoords":
        OPTIONS.bed = os.path.abspath(OPTIONS.bed)
        OPTIONS.chromosome_lengths = os.path.abspath(OPTIONS.chromosome_lengths)
        bed_parser = sv_to_bed.BedParser(OPTIONS.bed, OPTIONS.chromosome_lengths)
        bed_parser.read_bed()
        bed_parser.read_chromosomes()
        bed_parser.randomise_bed(int(OPTIONS.n_shuffles))

    elif OPTIONS.cmd == "getchromlens":
        OPTIONS.fasta = os.path.abspath(OPTIONS.fasta)
        chromosome_lengths_dict = sv_to_bed.get_chromosome_lengths(OPTIONS.fasta)

        for key, val in chromosome_lengths_dict.items():
            print(key+"\t"+str(val))

    elif OPTIONS.cmd == "plink2ldhat":
            OPTIONS.ped = os.path.abspath(OPTIONS.ped)
            OPTIONS.map = os.path.abspath(OPTIONS.map)

            if OPTIONS.fasta:
                OPTIONS.fasta = os.path.abspath(OPTIONS.fasta)
                chromosome_lengths_dict = sv_to_bed.get_chromosome_lengths(OPTIONS.fasta)
            else:
                OPTIONS.chromosome_lengths = os.path.abspath(OPTIONS.chromosome_lengths)
                chromosome_lengths_dict = \
                    sv_utilities.read_chromosome_lengths(OPTIONS.chromosome_lengths)

            plink_ped = sv_utilities.PlinkPed(OPTIONS.ped, OPTIONS.map,
                                              chromosome_lengths_dict)
            plink_ped.read_ped()
            plink_ped.read_map()
            plink_ped.filter_segregating()

            if OPTIONS.ldhot:
                plink_ped.to_ldhot()
                plink_ped.filter_missing()
                plink_ped.print_per_chromosome(ldhot=True)

            else:
                plink_ped.print_per_chromosome()

    elif OPTIONS.cmd == "watfsites":
        OPTIONS.sites = os.path.abspath(OPTIONS.sites)
        OPTIONS.locs = os.path.abspath(OPTIONS.locs)
        watf_sites = sv_utilities.finite_sites_watterson(OPTIONS.sites, OPTIONS.locs)

    else:
        OPTIONS.bubbles = os.path.abspath(OPTIONS.bubbles)
        bubbles = sv_gene_centric.PanGeneBubbles(OPTIONS.bubbles)
        bubbles.parse_bubbles()
        bubbles.get_missing_info()

    shutil.rmtree(OPTIONS.temp_dir)

main()
