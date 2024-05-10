#!/usr/bin/env python3

import sys
import os
import bdsg
from .sv_utilities import CactusGraph, verbose_function
from .sv_utilities import check_sv
import random

class VcfSvPositions:

    def __init__(
        
                 self,
                 vcf: str,
                 reference_name: str,
                 graph: str,
                 temp_dir: str,
                 vg_path=os.path.dirname(os.path.abspath(__file__))+"/../bin",
                 threads=1,
                 consider_others=False,
                 verbose=False,
                 sv_size=50

                 ):
        
        self.consider_others = consider_others
        self.vcf = vcf
        self.vcf_positions_dict = dict()
        self.reference_name = reference_name
        self.vg_path = vg_path
        self.threads = threads
        self.graph = CactusGraph(graph, temp_dir, vg_path = self.vg_path, threads = self.threads, verbose = verbose)
        self.graph.get_handlegraph()
        self.packed_graph = self.graph.packed_graph
        self.posov = self.graph.posov
        self.complex_variants = set()
        self.genotypes_dict = dict()
        self.verbose = verbose
        self.sv_size = sv_size

    # A function to parse the VCF file and get the positions of the SVs in the graph.
    def parse_vcf(self):

        # Print to stderr that the VCF file is being parsed.
        verbose_function(self.verbose, "Parsing VCF to get variant positions.")

        # A function to check if a variant is complex.
        def check_complex(genotypes, is_sv):

            if is_sv:

                for i in genotypes:

                    if i.isdigit():

                        if int(i) > 1:
                            return True
                        
            return False


        # A function to check if the genotype name is present in both the start and end positions.
        def check_both_present(genotype_name, start_positions, end_positions):

            if genotype_name in start_positions and genotype_name in end_positions:
                return True

            return False
        
        # A function to update the positions dictionary. To be used in a lambda function later.
        def update_positions(x, positions_dict, packed_graph, posov, start_handle, start):

            if start:
                path_name = packed_graph.get_path_name(packed_graph.get_path_handle_of_step(x)).split("#")[0]
                position = posov.get_position_of_step(x) + packed_graph.get_length(start_handle)
                positions_dict[path_name] = position

            else:
                path_name = packed_graph.get_path_name(packed_graph.get_path_handle_of_step(x)).split("#")[0]
                position = posov.get_position_of_step(x)
                positions_dict[path_name] = position

            return True

        # Open the VCF file.
        with open(self.vcf) as f:
            
            # For each line in the VCF file.
            for line in f:

                if line[0] == "#":

                    if line[0:6] == "#CHROM":
                        self.alt_names = line.strip().split("\t")[9:]
                
                # If the line is not a header line.
                else:
                    # Get the information from the line for that variant.
                    fields = line.strip().split("\t")
                    (start_id, end_id) = list(map(int, fields[2].split(">")[1:]))
                    chromosome = fields[0]
                    start_handle = self.graph.packed_graph.get_handle(start_id) # Use the gbz graph to get the start and end handles.
                    end_handle = self.graph.packed_graph.get_handle(end_id)
                    is_an_sv = check_sv(fields[3], fields[4].split(","), self.sv_size)
                    genotypes = fields[9:]
                    complex_variant = check_complex(genotypes, is_an_sv)

                    if chromosome not in self.genotypes_dict:
                        self.genotypes_dict[chromosome] = dict()

                    self.genotypes_dict[chromosome][fields[2]] = genotypes

                    if complex_variant:
                        self.complex_variants.add(fields[2])

                    start_positions = dict() # Empty dictionaries to store the positions of the SVs.
                    end_positions = dict()

                    # For each step on the start and end handles, update the positions dictionaries.                                                
                    self.packed_graph.for_each_step_on_handle(
                            
                                                        start_handle,
                                                           
                                                        lambda x:
                                                        update_positions(x, start_positions, self.packed_graph, self.posov, start_handle, True)
                                                           
                                                        )
                        
                    self.packed_graph.for_each_step_on_handle(
                            
                                                        end_handle,
                                                           
                                                        lambda x:
                                                        update_positions(x, end_positions, self.packed_graph, self.posov, start_handle, False)
                                                           
                                                        )
                                                
                        
                    for each in [self.reference_name] + self.alt_names:
                        
                        # Check if the genotype has both start and end positions.
                        if check_both_present(each, start_positions, end_positions):
                            
                            # If only svs are to be considered, check if the variant is an SV.
                            if not self.consider_others:

                                # If it is an SV, update the positions dictionary.
                                if is_an_sv:

                                    if each not in self.vcf_positions_dict:
                                        self.vcf_positions_dict[each] = dict()
                                    
                                    if chromosome not in self.vcf_positions_dict[each]:
                                        self.vcf_positions_dict[each][chromosome] = [[start_positions[each], end_positions[each], ">"+str(start_id)+">"+str(end_id)]]
                                    
                                    else:
                                        self.vcf_positions_dict[each][chromosome].append([start_positions[each], end_positions[each], ">"+str(start_id)+">"+str(end_id)])

                            else:
                                
                                # Otherwise, update the positions dictionary.
                                if each not in self.vcf_positions_dict:
                                    self.vcf_positions_dict[each] = dict()
                                    
                                if chromosome not in self.vcf_positions_dict[each]:
                                    self.vcf_positions_dict[each][chromosome] = [[start_positions[each], end_positions[each], ">"+str(start_id)+">"+str(end_id)]]
                                    
                                else:
                                    self.vcf_positions_dict[each][chromosome].append([start_positions[each], end_positions[each], ">"+str(start_id)+">"+str(end_id)])
                                    

    def print_bed(self):

        for each in [self.reference_name] + self.alt_names:
            
            for chromosome, positions in self.vcf_positions_dict[each].items():
                
                for positions_list in positions:
                    print(chromosome+"\t"+str(positions_list[0]-1)+"\t"+str(positions_list[1])+"\t"+positions_list[2]+";"+each)
                

class BedParser:

    def __init__(
        
                  self,
                  bed: str,
                  chromosome_lengths: str

                  ):
        self.bed = bed
        self.bed_dict = dict()
        self.chromosome_lengths = chromosome_lengths
        self.chromosome_lengths_dict = dict()
        
    def read_bed(self):
        
        with open(self.bed) as f:

            for line in f:
                fields = line.strip().split("\t")
                
                if fields[0] not in self.bed_dict:
                    self.bed_dict[fields[0]] = [int(fields[2]) - int(fields[1])]

                else:
                    self.bed_dict[fields[0]].append(int(fields[2]) - int(fields[1]))

    def read_chromosomes(self):

        with open(self.chromosome_lengths) as f:

            for line in f:
                fields = line.strip().split()
                self.chromosome_lengths_dict[fields[0]] = int(fields[1])

    def randomise_bed(
                       self,
                       n_shuffles: int

                      ):
        sys.stderr.write("Randomising BED file " + self.bed + ".\n")
        
        for i in range(n_shuffles):

            for chromosome, lengths in self.bed_dict.items():

                for length in lengths:
                    max_length = self.chromosome_lengths_dict[chromosome]
                    random_end = random.randint(0 + length, max_length)
                    random_start = random_end - length
                    print(chromosome + "\t" + str(random_start) + "\t" + str(random_end) + "\trandomisation_" + str(i))

def get_chromosome_lengths(genome_file: str):
    chromosome_lengths_dict = dict()
    
    with open(genome_file) as f:

        for line in f:

            if line[0] == ">":
                chromosome = line.strip().split()[0][1:]
                chromosome_lengths_dict[chromosome] = 0

            else:
                chromosome_lengths_dict[chromosome] += len(line.strip())

    return chromosome_lengths_dict
        
