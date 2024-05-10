#!/usr/bin/env python3

import re

class PanGeneBubbles:

    """

    A class for working with the output of pangene, which creates a file
    called 'bubbles.txt' using the k8 library.

    Attributes:
    bubbles: str
        The path to the bubbles.txt file
    bubbles_dict: dict
        A dictionary representation of the bubbles.txt file

    """

    def __init__(self,
                 bubbles: str):
        self.bubbles = bubbles
        self.bubbles_dict = dict()

    def parse_bubbles(self):

        """

        Parse the bubbles.txt file and store the information in a dictionary.

        """
        
        with open(self.bubbles, 'r') as f:

            for line in f:
                fields = line.strip().split("\t")

                if fields[0] == "CC":
                    continue

                elif fields[0] == "BB":
                    bubble_id = fields[1]

                    if bubble_id not in self.bubbles_dict:
                        self.bubbles_dict[bubble_id] = dict()

                    self.bubbles_dict[bubble_id]["start"] = fields[4]
                    self.bubbles_dict[bubble_id]["end"] = fields[5]
                    self.bubbles_dict[bubble_id]["ngenes"] = fields[7]
                    self.bubbles_dict[bubble_id]["nalleles"] = fields[6]
                    self.bubbles_dict[bubble_id]["genes"] = fields[8].split(",")

                elif fields[0] == "AL":

                    if "haplotype_walks" not in self.bubbles_dict[bubble_id]:
                        self.bubbles_dict[bubble_id]["haplotype_walks"] = []

                    if "haplotype_allele_counts" not in self.bubbles_dict[bubble_id]:
                        self.bubbles_dict[bubble_id]["haplotype_allele_counts"] = []

                    if "haplotype_individuals" not in self.bubbles_dict[bubble_id]:
                        self.bubbles_dict[bubble_id]["haplotype_individuals"] = []

                    self.bubbles_dict[bubble_id]["haplotype_walks"].append(fields[2])
                    self.bubbles_dict[bubble_id]["haplotype_allele_counts"].append(fields[1])
                    self.bubbles_dict[bubble_id]["haplotype_individuals"].append(fields[3])

                else:
                    continue

    def get_missing_info(self):

        """

        Based on the information in self.bubbles_dict, this will add information
        on how many consecutive genes are missing between alleles. For multiple
        alleles, this will be the maximum number of consecutive genes missing.

        Parameters:
        self

        Returns:
        None

        """
        def largest_consecutive_run(list1, list2):

            """

            Finds the largest consecutive run of elements in list 1 not in list 2 and vice versa.
            Returns the longest of the two.

            Parameters:
            list1: list
            list2: list

            Returns:
            int
            
            """

            list1_missing = []
            list2_missing = []

            for i in range(len(list1)):

                if list1[i] not in list2:
                    list1_missing.append(0)

                else:
                    list1_missing.append(1)

            for i in range(len(list2)):
                if list2[i] not in list1:
                    list2_missing.append(0)

                else:
                    list2_missing.append(1)

            max_consecutive = 0
            current_consecutive = 0

            for each in list1_missing:

                if each == 0:
                    current_consecutive += 1
                    max_consecutive = max(max_consecutive, current_consecutive)

                else:
                    current_consecutive = 0

            current_consecutive = 0

            for each in list2_missing:

                    if each == 0:
                        current_consecutive += 1
                        max_consecutive = max(max_consecutive, current_consecutive)
    
                    else:
                        current_consecutive = 0

            return max_consecutive    
        

        for bubble_id in self.bubbles_dict:
            haplotype_walks = self.bubbles_dict[bubble_id]["haplotype_walks"]
            haplotype_individuals = self.bubbles_dict[bubble_id]["haplotype_individuals"]
            max_consecutives = []
            max_consecutive_individuals = []

            for i in range(len(haplotype_walks)-1):

                for j in range(i+1, len(haplotype_walks)):
                    walk1 = re.split('<|>', haplotype_walks[i])
                    walk2 = re.split('<|>', haplotype_walks[j])
                    max_consecutive = largest_consecutive_run(walk1, walk2)
                    max_consecutives.append(max_consecutive)
                    max_consecutive_individuals.append((haplotype_individuals[i], haplotype_individuals[j]))

            max_consecutive_index = max_consecutives.index(max(max_consecutives))
            max_consecutive = max_consecutives[max_consecutive_index]
            max_consecutive_individuals = max_consecutive_individuals[max_consecutive_index]
            print(bubble_id+"\t"+self.bubbles_dict[bubble_id]["start"]+"\t"\
                  +self.bubbles_dict[bubble_id]["end"]+"\t"\
                  +self.bubbles_dict[bubble_id]["ngenes"]+"\t"
                  +self.bubbles_dict[bubble_id]["nalleles"]+"\t"\
                  +str(max_consecutive)+"\t"+max_consecutive_individuals[0]+"\t"\
                  +max_consecutive_individuals[1])
