---
title: "Find Paralogs"
author: "Ryan Peek"
date: "*Updated: 2018-08-01*"
output: 
  html_document:
    keep_md: true
    toc: yes
    toc_float: yes
    code_folding: hide
---



## Finding Paralogs with *ngsParalog*

ngsParalog is a program for detecting cryptic copy number variation from population-level, next generation sequencing (NGS) data. ngsParalog implements a likelihood method for estimating the probability of mismapping reads due to paralogs while modeling aspects of NGS data, which makes it effective at even very low sequencing depths. The program reads in pileup format data and outputs per site likelihood ratios of duplication. These likelihood ratios are asymptotically distributed as a 50/50 mixture of a Chi-Square(1 d.f.) and a point mass at zero. From this [page](https://github.com/tplinderoth/ngsParalog).

### Step 1: Check Quality and Make Bamlists

We can run quality filters/checks on our data and generate the preferred bamlists using other scripts/tools from prior parts of this pipeline.

### Step 2: Find SNPs 

Before moving forward we need to identify SNPs in our data. This script uses `angsd` to identify SNPs from a bamlist of samples, or a list of bamlists.

**01_get_snps.sh**

```bash

#!/bin/bash

#mkdir results_snp  ### All output goes here ###

infile=$1 ### list containing population names ###
ref='final_contigs_300.fasta' ### reference fasta file used in alignment ###
n=$(wc -l $infile | awk '{print $1}')

### Calculate maf and snps ###

x=1
while [ $x -le $n ] 
do

        pop=$(sed -n ${x}p $infile)

                echo "#!/bin/bash" > ${pop}.sh
                echo "" >> ${pop}.sh
                echo "" >> ${pop}.sh
                echo "#SBATCH -o /home/rapeek/projects/rasi_hybrid/slurm_outs/02_get_snps-%j.out" >> ${pop}.sh
                echo "" >> ${pop}.sh
                echo "angsd -bam bamlists/${pop}_25k_thresh.bamlist -out results_snp/$pop -ref $ref -GL 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 10 -minQ 20" >> ${pop}.sh
                echo "gunzip results_snp/${pop}*.gz" >> ${pop}.sh
                echo "cut -d$'\t' -f1-2  results_snp/${pop}.mafs | sed 1d > results_snp/${pop}.snp.pos" >> ${pop}.sh
                echo "echo 'finished running get_snps'" >> ${pop}.sh

                sbatch -J rp_snps -t 4800 --mem=60G -c 1 ${pop}.sh

        x=$(( $x + 1 ))

done
date

```


### Step 3: Find Paralogs

We need to trim/filter reads that align to more than one area in the genome (paralogs). Or identify through loci depth, where reads may be coming from two paralogs in in sequencing panel, but your reference genome assembled only one of those paralogous sequences. This would lead to a bias in the depth of the locus compared to the other (one has many more mapped than the other. This script does that using a likelihood ratio test.

**02_find_paralogs.sh**

```bash

#!/bin/bash

mkdir results_paralogs  ### All output goes here ###

infile=$1 ### list containing population names ###
ref='final_contigs_300.fasta' ### reference fasta file used in alignment ###
n=$(wc -l $infile | awk '{print $1}')

### Calculate paralog probabilities and get a list of paralogous loci  ###

x=1
while [ $x -le $n ] 
do

        pop=$(sed -n ${x}p $infile)

                echo "#!/bin/bash" > ${pop}_paralog.sh
                echo "" >> ${pop}_paralog.sh
                echo "#SBATCH -o out_slurms/05b_paralogs-%j.out" >> ${pop}_paralog.sh
                echo "" >> ${pop}_paralog.sh

                echo "samtools mpileup -b bamlists/${pop}_25k_thresh.bamlist -l results_snp/${pop}.snp.pos -f ${ref} > results_paralogs/${pop}.depth" >> ${pop}_paralog.sh
                echo "~/bin/ngsParalog/ngsParalog calcLR -infile results_paralogs/${pop}.depth > results_paralogs/${pop}.paralogs" >> ${pop}_paralog.sh
                echo "awk '(\$5 > '10')' results_paralogs/${pop}.paralogs | cut -c1-7 | uniq > results_paralogs/${pop}.paralogs.list" >> ${pop}_paralog.sh
                echo "grep '>' $ref | cut -c2- | grep -v -f results_paralogs/${pop}.paralogs.list | sed 's/$/:/' > results_paralogs${pop}.loci" >> ${pop}_paralog.sh
                sbatch -J rpparlog -t 4800 --mem=16G -p high -c 1 ${pop}_paralog.sh

        x=$(( $x + 1 ))

done

```


### Step 4: Do SFS (or something else)

Now we can create a site frequency spectrum output using our filtered data. Note, here we add a `-rf` flag which specifies the "good" loci we want to restrict our angsd analysis to, as well as the baits, to specify we want to use only filtered (no paralogs) loci from the sites list.

**03_makeSFS.sh**


```bash

#!/bin/bash -l

# mkdir results_sfs

infile=$1 ### list containing population names ###
filtered=$2 ## yes or no
ref='final_contigs_300.fasta' ### reference fasta file used in alignment, IMPORTANT!!! Reference must correspond to ancestral states, if not supply a different fasta file for -anc!!! ###

n=$(wc $infile | awk '{print $1}')

#<<'Comment' # for multiline comment

### Calculate saf files and the ML estimate of the sfs using the EM algorithm for each population ###
x=1
while [ $x -le $n ] 
do

        pop=$(sed -n ${x}p $infile)

                echo "#!/bin/bash -l" > ${pop}_sfs.sh
                echo "" >> ${pop}_sfs.sh
                echo "#SBATCH -o out_slurms/05a_unfolded_sfs-%j.out" >> ${pop}_sfs.sh
                echo "#SBATCH -J rpsfs" >> ${pop}_sfs.sh
                echo "" >> ${pop}_sfs.sh

                if [ "$2"="yes" ]; then
                        # if filter_paralogs has been run, use this:
                        echo "angsd -bam bamlists/${pop}_25k_thresh.bamlist -ref $ref -anc $ref -sites bait_lengths.txt -rf results_snp/${pop}.loci -out results_sfs/$pop -GL 2 -doSaf 1 -minMapQ 10 -minQ 20" >> ${pop}_sfs.sh
                else
                        # if no filter paralogs, use this:
                        echo "angsd -bam bamlists/${pop}_25k_thresh.bamlist -ref $ref -anc $ref -sites bait_lengths.txt -out results_sfs/$pop -GL 2 -doSaf 1 -minMapQ 10 -minQ 20" >> ${pop}_sfs.sh
                fi
                # then get sfs
                echo "realSFS results_sfs/${pop}.saf.idx -maxIter 100 > results_sfs/${pop}.sfs" >> ${pop}_sfs.sh
                echo "~/scripts/plotSFS.R results_sfs/${pop}.sfs" >> ${pop}_sfs.sh

                sbatch -t 3600 --mem=16G -c 1 ${pop}_sfs.sh

        x=$(( $x + 1 ))

done

```
