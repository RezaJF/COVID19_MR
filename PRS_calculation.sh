#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
# Number of cpu cores required
#$ -pe smp 1
# RAM requirement per cpu core
#$ -l h_vmem=32G

if [ $# -ne 3 ]
then
    echo "Usage: `basename $0` PATH_to_PRSice TARGET_BED OUTPUT_DIRECTORY"
    exit 1
fi

echo "`date -u` : Run started successfully" > logfile.txt


PATH_to_PRSice=$1
target=$2
out=$3


cat list_of_summary_statistic.txt | while read i
	do

		gunzip ${i}.gz
		awk -F "\t" 'BEGIN {OFS = ":"} {print $1,$2}' $i > xx

		paste $i xx > ${i}.assoc





		${PATH_to_PRSice} \
			--A1 ALT \
               		--A2 REF \
               		--bar-levels 1e-40,1e-30,1e-20,1e-10,1e-09,1e-08,1e-07,1e-06,1e-05,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.05,0.1,1 \
               		--base ${i}.assoc \
               		--beta  \
               		--binary-target F \
               		--bp BP \
               		--chr CHR \
               		--clump-kb 250 \
               		--clump-p 1.000000 \
               		--clump-r2 0.100000 \
               		--interval 5e-05 \
               		--lower 5e-08 \
               		--missing MEAN_IMPUTE \
               		--model add \
               		--out ${out}/${i}_prs \
               		--perm 10000 \
               		--pvalue P-value \
               		--score avg \
               		--seed 910185829 \
               		--snp CHR:BP \
               		--stat Effect \
               		--target ${target} \
               		--thread 1 \
               		--upper 0.5

rm ${i}.assoc xx
gzip ${i}


	done
	
mkdir PRS_scores
mv *_prs.* PRS_scores
