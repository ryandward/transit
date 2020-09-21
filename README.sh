#!/usr/bin/env bash

# for the first walkthrough, the chromosome we will use is the E. coli MG1655 chromosome
GENOME_TARGET="U00096.3"

conda create -n transit -c bioconda 'bowtie=1.2.3' entrez-direct git 'transit==3.1.0' 'samtools==1.10' 'bedtools==2.29.2' 'seqkit==0.11.0' 'exonerate==2.4.0' tabulate 'python>3'

conda activate transit

# fetch the nucleotide fasta from NCBI
esearch -db nuccore  -query $GENOME_TARGET |
  efetch -format fasta > ${GENOME_TARGET}.fasta &&
    file ${GENOME_TARGET}.fasta |
      grep -iq ascii &&
        echo "${GENOME_TARGET}.fasta contains data, continue to next step." ||
          echo "Emtpy file: ${GENOME_TARGET}.fasta, try efetch step again."

# fetch the feature table from NCBI
esearch -db nuccore  -query $GENOME_TARGET |
  efetch -format ft > ${GENOME_TARGET}.tbl &&
    file ${GENOME_TARGET}.tbl |
      grep -iq ascii &&
        echo "${GENOME_TARGET}.tbl contains data, continue to next step." ||
          echo "Emtpy file: ${GENOME_TARGET}.tbl, try efetch step again."

# build a bowtie index
bowtie-build -o1 ${GENOME_TARGET}.fasta $GENOME_TARGET

######################################################################
# now, find essential genes using transit

# use bowtie to align the reads to the genome
bowtie -S -m2 --threads 12 --no-unal --best --strata --tryhard \
$GENOME_TARGET <(awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 25) {print header, seq, qheader, qseq}}' ${GENOME_TARGET}_Reads/*.fastq) |
samtools view -S -b -F 4 > ${GENOME_TARGET}_aln.bam

# sort the bam then save it to a temp file
samtools sort ${GENOME_TARGET}_aln.bam > ${GENOME_TARGET}_aln.bam.sort

# move the temp file back to the original bam, to save space
mv ${GENOME_TARGET}_aln.bam.sort ${GENOME_TARGET}_aln.bam

# use bam2wig.pl hopefully saved in a perlbrew environment, more documentation later
bam2wig.pl  --in ${GENOME_TARGET}_aln.bam --start --rpm --out ${GENOME_TARGET}_aln.wig

# use this ad-hoc script to convert the feature table to the the custom format required by the authors of transit (.tbl -> .bed -> .prot_table)
./tbl2bed.awk ${GENOME_TARGET}.tbl |
awk '$7=="CDS" {print $1,$2,$3,$6,($3-$2)/3,$7,$8,$5,$4}' > ${GENOME_TARGET}.prot_table

# finally, run the algorithm to detect essential genes
transit tn5gaps -m2 -r mean ${GENOME_TARGET}_aln.wig ${GENOME_TARGET}.prot_table ${GENOME_TARGET}_transit_tn5gaps.tsv

##########################################################
# once two genomes have been compared, test for enrichment
# build fasta amino acid files to test for orthology

#! todo: add options to pile together multiple chromosomes, such as for plasmids

# make a directory to put all the intermediate files
mkdir Working

# first set two targets for comparison
GENOME_TARGET_1="U00096.3"
GENOME_TARGET_2="CP008706.1"

# convert the feature tables into bed files
./tbl2bed.awk ${GENOME_TARGET_1}.tbl |
awk '$7 == "CDS" && $8 == "complete"' > Working/${GENOME_TARGET_1}_coding.bed

./tbl2bed.awk ${GENOME_TARGET_2}.tbl |
awk '$7 == "CDS" && $8 == "complete"' > Working/${GENOME_TARGET_2}_coding.bed

# convert the bedfiles into fasta nucleotide sequences
bedtools getfasta -nameOnly -s -bed Working/${GENOME_TARGET_1}_coding.bed -fi ${GENOME_TARGET_1}.fasta |
sed 's/(.)//g'> Working/${GENOME_TARGET_1}_coding.fna

bedtools getfasta -nameOnly -s -bed Working/${GENOME_TARGET_2}_coding.bed -fi ${GENOME_TARGET_2}.fasta |
sed 's/(.)//g' > Working/${GENOME_TARGET_2}_coding.fna

# some nucleotide sequences are discontiguous, so merge them
awk 'FNR%2==1 && !($0 in seen){print; seen[$0]; getline; print} FNR%2==1 && $0 in seen {getline; print}' Working/${GENOME_TARGET_1}_coding.fna |
awk 'BEGIN{ORS=""} /^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' > Working/tmp && mv Working/tmp Working/${GENOME_TARGET_1}_coding.fna

awk 'FNR%2==1 && !($0 in seen){print; seen[$0]; getline; print} FNR%2==1 && $0 in seen {getline; print}' Working/${GENOME_TARGET_2}_coding.fna |
awk 'BEGIN{ORS=""} /^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' > Working/tmp && mv Working/tmp Working/${GENOME_TARGET_2}_coding.fna

# translate them into amino acid sequences
seqkit translate -f1 -x -w0 -T11 -x -M Working/${GENOME_TARGET_1}_coding.fna > Working/${GENOME_TARGET_1}_coding.faa

seqkit translate -f1 -x -w0 -T11 -x -M Working/${GENOME_TARGET_2}_coding.fna > Working/${GENOME_TARGET_2}_coding.faa

# use exonerate to generate a list of pair-wise best matches
# output the query, target, percent similarity, raw score
exonerate --bestn 1 Working/${GENOME_TARGET_1}_coding.faa Working/${GENOME_TARGET_2}_coding.faa --showalignment 0 --showvulgar 0 --ryo "%qi\t%ti\t%ps\t%s\t"\#ortho"\n" > Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv

exonerate --bestn 1 Working/${GENOME_TARGET_2}_coding.faa Working/${GENOME_TARGET_1}_coding.faa --showalignment 0 --showvulgar 0 --ryo "%ti\t%qi\t%ps\t%s\t"\#ortho"\n" >> Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv

# find orthos that occur more than once in the list, i.e. pair-wise best
awk '$5=="#ortho" {seen[$1 OFS $2 OFS $3 OFS $4]++} END{for(i in seen){if(seen[i]>1){print i}}}' Working/U00096.3_CP008706.1_orthos.tsv > Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best.tsv

# retrieve gene extended gene information
esearch -db gene -query "$GENOME_TARGET_1\[NUCL_ACCN\]" |
efetch -format tabular | awk 'FNR==1{print $3, $6, $7; seen[$3]} !($3 in seen) {seen[$3]; split($7, syn, ", "); for(i in syn) {print $3, $6, syn[i], $5;}}' > Working/${GENOME_TARGET_1}_genes.tsv

esearch -db gene -query "$GENOME_TARGET_2\[NUCL_ACCN\]" |
efetch -format tabular | awk 'FNR==1{print $3, $6, $7; seen[$3]} !($3 in seen) {seen[$3]; split($7, syn, ", "); for(i in syn) {print $3, $6, syn[i], $5;}}'> Working/${GENOME_TARGET_2}_genes.tsv

mkdir Enrichment 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos 2>/dev/null

# add gene IDs
awk 'NR == FNR {def[$3]=$1 OFS $2} NR != FNR && $1 in def{print $0, def[$1]}' Working/${GENOME_TARGET_1}_genes.tsv Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best.tsv  | sort -k6,6V > ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv

# only use gene IDs
awk '{print $6}' Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv > Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv

# number of genes
BACKGROUND="$(wc -l Working/${GENOME_TARGET_1}_coding.bed | cut -d' ' -f1)"

# find enrichment of orthologs
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF -m=normal -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv -t=sym -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos/ -db=GO_all -p=yes -bg=$BACKGROUND -org=ecocyc

# compare essentiality of both organisms wrt the first
awk 'NR==FNR {ess[$1]=$NF} NR != FNR {print $0, ess[$1]}' ${GENOME_TARGET_1}_transit_tn5gaps.tsv Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv | awk 'NR==FNR {ess[$1]=$NF} NR != FNR {print $0, ess[$2]}' ${GENOME_TARGET_2}_transit_tn5gaps.tsv - | awk '{print $1, $7, $2, $8, $5, $6}'> Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv

# create empty directories
mkdir Enrichment 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN 2>/dev/null

# find false positives genome 2 wrt 1
awk '$2=="Non-essential" && $4=="Essential" {print $NF}' Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv > Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP.tsv

# find false negatives genome 2 wrt 1
awk '$2=="Essential" && $4=="Non-essential" {print $NF}' Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv > Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN.tsv

# number of orthologs
BACKGROUND="$(wc -l Working/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv | cut -d' ' -f1)"

# find false negative enrichment
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF -m=normal -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN.tsv -t=sym -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN/ -db=GO_all -p=yes -bg=$BACKGROUND -org=ecocyc

# find false positive enrichment
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF -m=normal -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP.tsv -t=sym -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP/ -db=GO_all -p=yes -bg=$BACKGROUND -org=ecocyc

conda deactivate
