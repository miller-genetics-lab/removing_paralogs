#!/bin/bash

mkdir results_snp  ### All output goes here ###

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
                echo "#SBATCH -o /home/slurm_outs/02_get_snps-%j.out" >> ${pop}.sh
                echo "" >> ${pop}.sh
                echo "angsd -bam bamlists/${pop}_25k_thresh.bamlist -out results_snp/$pop -ref $ref -GL 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 10 -minQ 20" >> ${pop}.sh
                echo "gunzip results_snp/${pop}*.gz" >> ${pop}.sh
                echo "cut -d$'\t' -f1-2  results_snp/${pop}.mafs | sed 1d > results_snp/${pop}.snp.pos" >> ${pop}.sh
                echo "echo 'finished running get_snps'" >> ${pop}.sh

                sbatch -J rp_snps -t 4800 --mem=60G -c 1 ${pop}.sh

        x=$(( $x + 1 ))

done
date