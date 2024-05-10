#!/usr/bin/env python3

import math
import re
import textwrap
from typing import List, Dict, Tuple
import subprocess
import numpy
import os
import sys
import bisect
import bdsg # Requires bdsg to be installed and available in the Python path.
import concurrent.futures
import itertools

class CactusGraph:

    def __init__(

                  self,
                  graph: str,
                  temp_dir: str,
                  graph_type = "gbz", 
                  threads = 1, 
                  vg_path = os.path.dirname(os.path.abspath(__file__))+"/../bin",
                  verbose = False
                  
                  ):
        
        self.graph = graph 
        self.packed_graph = bdsg.bdsg.PackedGraph() 
        self.posov = bdsg.bdsg.PackedPositionOverlay() 
        self.graph_type = graph_type 
        self.threads = str(threads) 
        self.vg_path = vg_path + "/vg" 
        self.vcf_file = "temp.vcf"
        self.temp_dir = temp_dir
        self.verbose = verbose

    def get_handlegraph(self): 
        verbose_function(self.verbose, "Making handlegraph with " + self.graph + ".")
        
        if self.graph_type == "gbz": 

            with open(self.temp_dir+"/"+"temp.pg", "wb") as temp_file: 
                packed_graph = subprocess.Popen([
                    
                                                  self.vg_path, 
                                                  "convert",
                                                  "--packed-out", 
                                                  "--threads",
                                                  self.threads,
                                                  self.graph
                                                  
                                                  ],

                                                stdout = temp_file)
                packed_graph.wait()
                self.packed_graph.deserialize(self.temp_dir+"/"+"temp.pg") 
                self.posov = bdsg.bdsg.PackedPositionOverlay(self.packed_graph) 
                self.graph = self.temp_dir+"/"+"temp.pg" 

        elif self.graph_type == "pg": 
            self.packed_graph.deserialize(self.graph) 
            self.posov = bdsg.bdsg.PackedPositionOverlay(self.packed_graph)

        else:
            sys.stderr.write("Graph type is one of either 'pg' or 'gbz'. Current graph type is '" + self.graph_type + "'\n.")
            sys.exit()
            
    def get_vcf(self, ref_path: str): 
        sys.stderr.write("Making VCF from " + self.graph + ".\n")
        
        if self.graph_type != "gbz": 
            sys.stderr.write("Need gbz output from cactus-pangenome to create vcf. Graph type is '" + self.graph_type + "'.\n")
            sys.exit()

        else: 
                
            with open(self.temp_dir+"/"+self.vcf_file, "w") as temp_file:
                deconstruct_gbz = subprocess.Popen([
                        
                                                    self.vg_path,
                                                    "deconstruct",
                                                    "--threads",
                                                    self.threads,
                                                    "--path-traversals",
                                                    "--ploidy",
                                                    "1",
                                                    "--path-prefix",
                                                    ref_path,
                                                    self.graph

                                                    ],

                                                    stdout = temp_file ,
                                                    text=True )
                deconstruct_gbz.wait()
    
    def get_graph_paths(self):        
        sys.stderr.write("Getting all path names from graph.\n")
        all_path_names = set()
        
        self.posov.for_each_path_handle(lambda x:
                                                all_path_names.add(self.posov.get_path_name(x)) or True)

        self.all_path_names = dict()
        
        for path_name in all_path_names:
            genotype = path_name.split("#")[0]
            chromosome = path_name.split("#")[2].split("[")[0]
            
            if chromosome not in self.all_path_names:
                self.all_path_names[chromosome] = dict()

            if genotype not in self.all_path_names[chromosome]:
                self.all_path_names[chromosome][genotype] = set()

            self.all_path_names[chromosome][genotype].add(path_name)
        

class CactusVcf: 

    def __init__(

                 self,
                 all_path_names: dict,
                 vcf_file: str, 
                 graph: bdsg.bdsg.PackedPositionOverlay, 
                 ref_path: str,
                 temp_dir: str,
                 packed_graph_file = "temp.pg",
                 threads = 1,
                 vg_path = os.path.dirname(os.path.abspath(__file__))+"/../bin", 
                 max_threads = 5,

                 ):
    
        self.vcf_file = vcf_file 
        self.graph = graph
        self.ref_path = ref_path
        self.complex_snarls = dict()
        self.packed_graph_file = packed_graph_file
        self.threads = str(threads)
        self.vg_path = vg_path + "/vg"
        self.max_threads = max_threads
        self.temp_dir = temp_dir
        self.all_path_names = all_path_names
        
    def get_snarl_info(self): 
        sys.stderr.write("Getting snarl info from " + self.vcf_file + ".\n")

        def check_alt_size(alts): 
            
            for alt in alts.split(","):
                
                if len(alt) >= 50:
                    return True
                
                return False
            
        with open(self.temp_dir+"/"+self.vcf_file) as vcf_file_handle: 
            
            for vcf_line in vcf_file_handle: 
                
                if not vcf_line[0] == "#": 
                    vcf_fields = vcf_line.strip().split("\t") 
                    chromosome = self.get_chromosome(vcf_fields[0]) 
                    variant_id = vcf_fields[2] 
                    snarl_id = tuple(map(int, variant_id.split(">")[1:3])) 
                    is_sv = False 

                    if len(vcf_fields[3]) < 50: 
                        is_sv = check_alt_size(vcf_fields[4])

                    else:
                        is_sv = True 

                    if is_sv: 
                        vcf_genotypes = vcf_fields[9:]
                                                
                        if max([int(x) for x in vcf_genotypes if x.isdigit()]) > 1:
                                    
                            if chromosome not in self.complex_snarls:
                                self.complex_snarls[chromosome] = [snarl_id] 
                            
                            else:
                                self.complex_snarls[chromosome].append(snarl_id) 
                else: 
                    
                    if vcf_line[0:6] == "#CHROM":
                        self.all_genotypes = [self.ref_path] + vcf_line.strip().split("\t")[9:]

    def get_chromosome(self, chromosome_field):

        if "#" in chromosome_field:
            return chromosome_field.split("#")[2]
            
        else:
            return chromosome_field
            
    def get_all_paths(self): 
        self.complex_snarl_positions = dict()

        for chromosome in self.complex_snarls.keys():
            
            self.complex_snarl_positions[chromosome] = dict()
            
            for genotype_name in self.all_genotypes:
                self.complex_snarl_positions[chromosome][genotype_name] = []

    
    def get_positions(self):
        sys.stderr.write("Getting snarl positions in all samples.\n")        
        
        self.complex_snarl_average_lengths = dict() 

        for chromosome, snarls in self.complex_snarls.items(): 

            for snarl in snarls:

                start_handle = self.graph.get_handle(snarl[0])
                end_handle = self.graph.get_handle(snarl[1])

                start_genotypes = set() 
                end_genotypes = set() 

                self.graph.for_each_step_on_handle(start_handle,
                                                       lambda x:
                                                           start_genotypes.add(
                                                                   self.graph.get_path_name(self.graph.get_path_handle_of_step(x)).split("#")[0]
                                                                   ) or True)

                self.graph.for_each_step_on_handle(end_handle,
                                                       lambda x:
                                                           end_genotypes.add(
                                                               self.graph.get_path_name(self.graph.get_path_handle_of_step(x)).split("#")[0]
                                                               ) or True)

                full_snarl_genotypes = start_genotypes & end_genotypes

                if (len(start_genotypes) != len(set(start_genotypes))):
                    sys.stderr.write("Duplicate on start.")
                    sys.stderr.write("Snarl ID: "+"-".join(map(str, list(snarl)))+"\n")
                    sys.exit()

                if (len(end_genotypes) != len(set(end_genotypes))):
                    sys.stderr.write("Duplicate on end.")
                    sys.stderr.write("Snarl ID: "+"-".join(map(str, list(snarl)))+"\n")
                    sys.exit()

                
                self.graph.for_each_step_on_handle( start_handle,
                                                    lambda x:
                                                            self.complex_snarl_positions[chromosome][   
                                                                                        self.graph.get_path_name(
                                                                                            self.graph.get_path_handle_of_step(x)
                                                                                            ).split("#")[0]
                                                                                        ].append(
                                                                                            self.graph.get_position_of_step(x) + self.graph.get_length(self.graph.get_handle_of_step(x))
                                                                                            ) or True )
                 
                self.graph.for_each_step_on_handle( end_handle,
                                                    lambda x:
                                                            self.complex_snarl_positions[chromosome][   
                                                                                        self.graph.get_path_name(
                                                                                            self.graph.get_path_handle_of_step(x)
                                                                                            ).split("#")[0]
                                                                                        ].append(
                                                                                            self.graph.get_position_of_step(x)+1
                                                                                            ) or True )
               
                total_lengths = []
                
                for genotype in full_snarl_genotypes:
                        total_lengths.append(self.complex_snarl_positions[chromosome][genotype][-1] - self.complex_snarl_positions[chromosome][genotype][-2])
                            
                self.complex_snarl_average_lengths[snarl] = numpy.sum(total_lengths) / len(total_lengths)


    def deconstruct_all_v_all(self, chromosome):
        
        all_paths = list(self.complex_snarl_positions[chromosome].keys())
        all_v_all_sv_counts = dict()  
        pairs = itertools.combinations(all_paths, 2)
        
        for pair in pairs:
            sys.stderr.write("Doing " + chromosome + " " + pair[0] + " " + pair[1] +".\n")

            with open(self.temp_dir+"/"+chromosome+"_temp_mod.pg", "wb") as temp_file:

                vg_command = [self.vg_path, "mod"]

                for path in self.all_path_names[chromosome][pair[0]]:
                    vg_command.extend(["--keep-path", path])

                for path in self.all_path_names[chromosome][pair[1]]:
                    vg_command.extend(["--keep-path", path])
                
                vg_command.append(self.temp_dir+"/"+self.packed_graph_file)
                vcf_combination = subprocess.Popen(vg_command, stdout = temp_file )
                vcf_combination.wait()

            dec_combination = subprocess.Popen([ self.vg_path,
                                                 "deconstruct",
                                                 "--threads",
                                                 self.threads,
                                                 self.temp_dir+"/"+chromosome+"_temp_mod.pg" ], stdout = subprocess.PIPE, text=True )

            vcf_iter = dec_combination.stdout.readlines()

            for vcf_line in vcf_iter:

                if vcf_line[0] != "#":
                    vcf_fields = vcf_line.strip().split("\t")
                    path_name = vcf_fields[0]
                    ref_position = int(vcf_fields[1])
                    snarl_positions = self.complex_snarl_positions[chromosome][path_name.split("#")[0]]
                    insertion_index = bisect.bisect_left(snarl_positions, ref_position)

                    if insertion_index % 2 != 0:
                        node_id = self.complex_snarls[chromosome][int((insertion_index-1)/2)]

                        if len(vcf_fields[3]) >= 50 or len(vcf_fields[4]) >= 50:

                            if node_id not in all_v_all_sv_counts:
                                all_v_all_sv_counts[node_id] = 1

                            else:
                                all_v_all_sv_counts[node_id] += 1

        return all_v_all_sv_counts
        
    def parallel_deconstruct_all_v_all(self):
        self.all_v_all_sv_counts = dict()
        n_processes = int(self.max_threads) // int(self.threads)
        chromosomes = list(self.complex_snarl_positions.keys())
    
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_processes) as executor:
            future_to_index = {executor.submit(self.deconstruct_all_v_all, element): index for index, element in enumerate(chromosomes)}

            for future in concurrent.futures.as_completed(future_to_index):
                index = future_to_index[future]
                result = future.result()
                self.all_v_all_sv_counts[chromosomes[index]] = result

    def get_chromosome_lengths(self, chromosome_lengths_file):
        self.chromosome_lengths = dict()

        with open(chromosome_lengths_file) as f:
            
            for line in f:
                self.chromosome_lengths[line.split()[0]] = line.strip().split()[1]

def check_sv(ref_allele, alt_alleles, sv_size=50):

    if len(ref_allele) >= sv_size:
        return True
        
    else:
            
        for i in alt_alleles:
                
            if len(i) >= sv_size:
                return True
                
    return False

def verbose_function(verbose, message):

    if verbose:
        sys.stderr.write(message + "\n")

def apply_jukes_cantor(n_sv_window):

    return -(4/3) * math.log(1-(3/4)*n_sv_window)

def reformat_vcf_egglib(vcf_file, temp_dir):

    """

    This function reformats a VCF file to be compatible with egglib. It is a bit of a workaround, as egglib does not support sites with over 10 alleles.

    vcf_file: The name of the VCF file to be reformatted.
    temp_dir: The temporary directory to write the reformatted VCF file to.

    Returns:

    None

    """
            
    with open(vcf_file) as f:

        with open(temp_dir+"/temp_egg.vcf", 'w') as of:
                
            for line in f:

                if line[0] != "#":
                    fields = line.strip().split()
                    fields[9:] = ["1" for i in fields[9:]]
                    of.write("\t".join(fields)+"\n")

                else:
                    of.write(line)

def split_vcf(vcf: str, chromosome: str, temp_dir: str):

    """

    This function splits a VCF file by chromosome.

    Parameters:
    
    vcf:        The name of the VCF file to be split.
    chromosome: The chromosome to split the VCF file by.
    temp_dir:   The temporary directory to write the split VCF file to.

    Returns:

    None

    """

    with open(temp_dir + "/" + chromosome + ".vcf", "w") as out_file:

        with open(vcf) as vcf_file:
            
            for line in vcf_file:

                if line.startswith("#"):
                    
                    if line.startswith("##contig"):

                        if line.split("=")[1].split(",")[0] == chromosome:
                            out_file.write(line)

                    else:
                        out_file.write(line)                              
            
                else:

                    if line.split()[0] == chromosome:
                        out_file.write(line)

def read_chromosome_lengths(chromosome_lengths_file: str) -> Dict[str, int]:

    """

    This function reads a file with chromosome lengths in it.

    Parameters:

    chromosome_lengths_file: The name of the file with chromosome lengths in it.

    Returns:

    chromosome_lengths:      A dictionary of chromosome names and their lengths.

    """

    chromosome_lengths = dict()

    with open(chromosome_lengths_file) as f:
        
        for line in f:
            chromosome_lengths[line.split()[0]] = int(line.strip().split()[1])

    return chromosome_lengths

class PlinkPed:

    """

    A class to hold information from a PLINK PED file.

    """

    def __init__(self, 
                 ped: str,
                 map: str,
                 chromosome_lengths: Dict[str, int]):
        self.ped = ped
        self.map = map
        self.seqs = dict()
        self.chromosome_lengths = chromosome_lengths
        self.filtered_sites = dict()
        self.map_numbers = dict()
        self.map_positions = []
        self.ldhot_seqs = dict()
        self.filtered_sites_missing = dict()

    def read_ped(self):
            
        """
    
        This function reads a PED file
        and creates a member variable called
        'seqs', which is a dictionary of sequence
        IDs and their alleles.

        """
        
        with open(self.ped) as f:

            # For each line the PED file, add
            # to seqs, seqs => seq_id => alleles.
            for line in f:
                fields = line.strip().split()
                seq_id = fields[0]
                self.seqs[seq_id] = fields[6::2]

    def read_map(self):

        """

        This function reads a MAP file
        and creates two member variables called
        'map_numbers' and 'map_positions', which
        are dictionaries of chromosome names and
        their corresponding map numbers and positions.

        """

        i = 0

        with open(self.map) as f:

                # For each line in the MAP file.
                for line in f:
                    fields = line.strip().split()
                    chromosome = fields[0] # First field is chromosome name.

                    # If the chromosome is not in the map_numbers dictionary.
                    if chromosome not in self.map_numbers:
                        self.map_numbers[chromosome] = [i] # Number is line no.

                    else:
                        self.map_numbers[chromosome].append(i)

                    i += 1 # Increment line number.

                    # map_positions is just a list of positions from start to end.
                    self.map_positions.append(fields[3])

    def to_ldhot(self):

        """

        A function to convert the sequences in the PED file
        to ldhot format (01).

        """

        # Dictionary to contain the alleles at each site.
        # Of the form site => alleles set.
        alleles_at_sites = dict()

        # Dictionary to contain a dictionary of alleles for each site.
        # Of the form site => allele => allele number.
        alleles_dicts_at_sites = dict()

        # For each sequence in the PED file.
        for key, val in self.seqs.items():

            # For each site in the sequence.
            for i in range(len(val)):

                # If the site is not in the alleles_at_sites dictionary.
                if i not in alleles_at_sites:
                    alleles_at_sites[i] = set() # Create an empty set for alleles.

                alleles_at_sites[i].add(val[i]) # Otherwise, add the allele to the set.

        # For each of the alleles_at_sites dictionary items.
        for key, val in alleles_at_sites.items():
            # Get a list of alleles that are not missing data.
            alleles_list = [i for i in list(val) if i in ["A", "C", "G", "T", "a", "c", "g", "t"]]

            # If this site is not in the alleles_dicts_at_sites dictionary.
            if key not in alleles_dicts_at_sites:
                alleles_dicts_at_sites[key] = dict() # Make an empty dictionary for alleles.

            # For all alleles at this site, in alphabetical order.
            for i in range(len(alleles_list)): # Allele points to list index, 0 or 1.
                alleles_dicts_at_sites[key][alleles_list[i]] = str(i)

        # Convert sequences to ldhot sequences.
        for key, val in self.seqs.items():
            self.ldhot_seqs[key] = []

            for i in range(len(val)):
                if val[i] in ["A", "C", "G", "T", "a", "c", "g", "t"]:
                    self.ldhot_seqs[key].append(str(alleles_dicts_at_sites[i][val[i]]))

                else:
                    self.ldhot_seqs[key].append("N")

    def filter_segregating(self):

        """

        This function filters out sites from the
        PED that are not segregating.

        """

        def check_segregating(map_numbers, seqs):
            seg_sites = []

            for i in map_numbers:
                alleles = set()

                for seq in seqs.values():

                    if seq[i] in ["A", "C", "G", "T", "a", "c", "g", "t"]:
                        alleles.add(seq[i])

                if len(alleles) > 1:
                    seg_sites.append(i)

            return seg_sites

        
        for key, val in self.map_numbers.items():
            seg_sites = check_segregating(val, self.seqs)
            self.filtered_sites[key] = seg_sites

    def filter_missing(self):

        """

        This function refines the filtered_sites dictionary and creates
        a new dictionary called filtered_sites_missing, which contains
        only segregating sites with no missing data.

        """

        def check_missing(map_numbers, seqs):
            non_missing = []

            for i in map_numbers:
                missing = False

                for seq in seqs.values():

                    if seq[i] not in ["A", "C", "G", "T", "a", "c", "g", "t"]:
                        missing = True

                if not missing:
                    non_missing.append(i)

            return non_missing

        for key, val in self.filtered_sites.items():
            non_missing = check_missing(val, self.seqs)
            self.filtered_sites_missing[key] = non_missing

    def print_per_chromosome(self, ldhot=False):

        """
        
        This function prints the alleles for each
        chromosome to a file.

        """

        def natural_sort_key(s):

            """
            
            This is just used to get sorted chromosome IDs based on chromosome number.

            """

            return [int(text) if text.isdigit() else text.lower()
                    for text in re.split(r'(\d+)', s)]


        if not ldhot:

            # self.filtered_sites is chromosome => filtered sites.
            for key, val in self.filtered_sites.items():
                # map_positions_chr is a list of filtered map positions.
                map_positions_chr = [self.map_positions[i] for i in val]

                with open(key+"_sites.txt", "w") as f: # Open a sites file for this chromosome.
                    # Write the number of sequences, number of filtered sites.
                    f.write(str(len(self.seqs))+" "+str(len(map_positions_chr))+" 1\n")

                # For each sequence in the PED file.
                for ikey, ival in self.seqs.items():
                    # seqs_chr_indv is a list of alleles for this sequence at the filtered sites.
                    seqs_chr_indv = [ival[i] for i in val]

                    # Append to the sites file for this chromosome.
                    with open(key+"_sites.txt", "a") as f:
                        f.write(">"+ikey[0:30]+"\n"+textwrap.fill("".join(seqs_chr_indv).replace("0","N"), width=1000)+"\n")

                # Open a locs file for this chromosome.
                with open(key+"_locs.txt", "w") as f:
                    chromosome_lengths_sorted = dict(sorted(self.chromosome_lengths.items(), 
                                                            key=lambda item: natural_sort_key(item[0])))
                    chromosome_lengths_sorted_print = dict()
                    i = 1

                    for chromosome in chromosome_lengths_sorted:
                        chromosome_lengths_sorted_print[str(i)] = self.chromosome_lengths[chromosome]
                        i += 1

                    f.write(str(len(map_positions_chr))+" "+str(float(chromosome_lengths_sorted_print[key])/float(1000))+" L\n")
                    f.write("\n".join([str(float(i)/float(1000)) for i in map_positions_chr])+"\n")

        else:

            for key, val in self.filtered_sites_missing.items():
                map_positions_chr = [self.map_positions[i] for i in val]

                with open(key+"_sites.txt", "w") as f:
                    f.write(str(len(self.ldhot_seqs))+" "+str(len(map_positions_chr))+" 1\n")
                for ikey, ival in self.ldhot_seqs.items():
                    seqs_chr_indv = [ival[i] for i in val]

                    with open(key+"_sites.txt", "a") as f:
                        f.write(">"+ikey[0:30]+"\n"+textwrap.fill("".join(seqs_chr_indv), width=1000)+"\n")

                with open(key+"_locs.txt", "w") as f:
                    chromosome_lengths_sorted = dict(sorted(self.chromosome_lengths.items(), 
                                                            key=lambda item: natural_sort_key(item[0])))

                    chromosome_lengths_sorted_print = dict()
                    i = 1

                    for chromosome in chromosome_lengths_sorted:
                        chromosome_lengths_sorted_print[str(i)] = self.chromosome_lengths[chromosome]
                        i += 1

                    f.write(str(len(map_positions_chr))+" "+str(float(chromosome_lengths_sorted_print[key])/float(1000))+" L\n")
                    f.write("\n".join([str(float(i)/float(1000)) for i in map_positions_chr])+"\n")

def finite_sites_watterson(sites, locs):

    """

    This function calculates a finite sites version of
    Watterson's theta using the sites and locs files used
    with LDhat.

    """

    def calculate_harmonic(n):

        """
        
        This function calculates the first term of
        the equation for the finite sites Watterson's
        theta.
    
        """

        harmonic = 0

        for i in range(1, n-1):
            harmonic += 1/i

        return harmonic ** -1

    with open(sites) as f:
        (n, S, undef) = f.readline().split()

    with open(locs) as f:
        (undef, L, undef) = f.readline().split()
    
    S = float(S)
    n = int(n)
    L = float(L)*1000

    print(calculate_harmonic(n) * math.log(L/(L-S)))
