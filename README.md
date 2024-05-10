# svstats: a package for analysing pan-genome graphs

svstats is a package for analysing pan-genome graphs. It relies on the bdsg API and vg software. So far, it has one main function, which is to calculate something called 'svpi' in sliding windows along a genome. This repository is a work in progress and will be regularly updated. There will eventually be more detailed installation instructions and instructions on running the program.

The package also has a few utilities for things. For example, one can create a TSV with the lengths of all chromosomes from a fasta with the sub-command 'getchromlens' or get random coordinates from a BED file with 'randcoords'. The sub-command 'svcoords' will give you the coordinates of all SVs in a VCF for each genotype, using a matching GBZ. The GBZ must have position information stored for all individuals, and individual names must match the VCF.

WARNING: I am not a professional software developer and the code here is my best attempt at solving problems related to my research. I have put the code in this repository to aid reproducibility and transparency but also because some of the functions may be useful to the pan-genomics research community. I am happy to address installation issues and bugs. While the software is broadly functional, it may not work with certain types of data (see input data descriptions below).

## Package overview

As I was analysing a pan-genome dataset, I came across various problems that I had to solve using scripting. The first of these problems was calculation of a new statistic called 'svpi'. This is a statistic inspired by the synteny diversity statistic of [Jiao et al. (2020)](https://doi.org/10.1038/s41467-020-14779-y), but works on counts of >= 50 bp variants between pairs of individuals, rather than the fraction of the region that is rearranged. In plain English, it is the average number of >= 50 bp variants between all individuals in a particular sliding window. The statistic is normalised to the average size of the sliding window across individuals. The window is set relative to the reference in the VCF, but for each individual, it will calculate the actual number of base pairs present in the window and take a window length average across individuals for normalising the svpi statistic.

After developing the scripts necessary for this, I came across several other problems that I also integrated into the package. Some are simple format conversion issues and others are also related to structural variant analysis in pan-genomes. In the following section, I will detail usage of the main command and the different sub-commands.

## Installation

The package requires vg 1.54.0 (https://github.com/vgteam/vg), egglib (https://egglib.org/) and the bdsg library (https://github.com/vgteam/libbdsg). The vg binary is provided with the svstats package in `bin`. Installation of bdsg is left to the user. Before running svstats, make sure bdsg is in PYTHONPATH:
```bash
export PYTHONPATH=/path/to/bdsg/lib:$PYTHONPATH
```

## The svstats packag sub-commands

The svstats package contains six sub-commands.

    svpi          Calculates svpi in user-defined sliding windows across a VCF from a minigraph cactus alignment.

    svcoords      Creates a BED file for each individual in the VCF with coordinates of all structural variants.

    randcoords    Creates BED files with random sets of coordinates for each individual. The random sets will
                  match the number and sizes of the structural variant coordinates.

    getchromlens  Creates a TSV file with chromosome lengths for all chromosomes in a FASTA file.

    plink2ldhat   Reformats plink-style .ped and .map formatted files to ldhat sites and locs files. If the switch
                  `--ldhot` is used, the format will be coded with zeros and ones and missing data will be removed.

    watfsites     Calculates the finite sites version of Watterson's theta using the formula from Auton et al.


### svpi

svpi requires the following things:

`reference_name`

This is the name of the reference sequence used in the vcf file provided. The script uses the sample names from the VCF and needs to add the reference name to this list of samples. The reference name must match one of the paths in the gbz.

`gbz`

This is a gbz file containing the genome graph used for finding the positions and lengths of variants relative to the reference genome. All individuals in this graph must be promoted to reference sense paths. The graph I used in my research was constructed by cactus minigraph with the following command:

```bash

references=($(cut -f 1 seqfile.txt | perl -pe "s/\n/ /g"))

cactus-pangenome ./jobStore \
                 seqfile.txt \
                 --outDir my_pg_dir \
                 --outName my_pg \
                 --reference ${references[@]} \
		        --gbz full clip filter 

```
All individuals in this graph are reference paths.
`vcf`
This is the gbz graph called against the named reference.
You can use vg to create this.
```bash
vg deconstruct --path-prefix my_ref \
               --path-traversals \
	           --ploidy 1 \
	           --contig-only-ref \
	           my_pg_dir/my_pg.full.gbz > my_pg_dir/my_pg.full.vcf

```

Unfortunately, the `svpi` statistic can only be calculated for haploid genomes or haplotypes. The genotypes in the VCF must be either '0', '1' or '.', and cannot be diploid calls, such as '1/1', '1|1', '0/0', '0|0', './.' or '.|.'. This may be fixed in a future version.

`vcf_files`

This is a directory containing VCF files for all possible pairwise comparisons between individuals in the pan-genome graph. The number of comparisons grows exponentially with the number of individuals, and starts to get very large beyond 30 genomes. One could in theory do all of these comparisons on a cluster. In my analyses, I had 25 genomes, and ran an all vs all comparison with cactus using the following command, assuming all assemblies are in a directory called `assemblies`:

```bash
mkdir seqfiles
assemblies_dir=`pwd`/assemblies
assemblies=(`ls assemblies`)

for ((i=0;i<${#assemblies[@]}-1;i++)); do

    for ((n=i+1;n<${#assemblies[@]};n++)); do

        file1=${assemblies[i]}
	    file2=${assemblies[n]}
	    strain1=${file1%.*}
	    strain2=${file2%.*}
	    seqfile=seqfiles/${i}_${n}_seqfile.txt
	
	    echo $strain1 $assemblies_dir/$file1 | perl -pe "s/ /\t/g" > $seqfile
	    echo $strain2 $assemblies_dir/$file2 | perl -pe "s/ /\t/g" >> $seqfile

	done
done

for seqfile in seqfiles/*; do

    prefix=${seqfile%.*}
    prefix=${prefix##*/}
    jobstore=${prefix}_js
    outdir=${prefix}_od
    outname=`cut -f 1 $seqfile  | perl -pe "s/\n/ /g" | awk '{print $1"_"$2}'`
    references=($(cut -f 1 $seqfile))

    cactus-pangenome $jobstore \
                 $seqfile \
                 --outDir $outdir \
                 --outName $outname \
                 --reference ${references[@]} \
		         --gbz full clip filter \
		         --vcf full clip filter
done
```

The names of the assembly files must match the names of the assemblies in `seqfile.txt`, used to create the master gbz and vcf (first bash command in this section). For example, the `seqfile.txt` used to create the master graph could have:
`strain1    /path/to/strain1.fasta`
Then, the same assembly file `strain1.fasta` could be put into the `assemblies` directory for the all against all comparison.
