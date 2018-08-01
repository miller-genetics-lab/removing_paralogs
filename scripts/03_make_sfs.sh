#!/bin/bash -l

mkdir results_sfs

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
