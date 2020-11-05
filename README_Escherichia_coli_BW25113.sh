#!/bin/bash
ORGANISM="Escherichia_coli_BW25113"
GENOME_TARGET=(CP009273.1) #https://www.ncbi.nlm.nih.gov/nuccore/CP009273.1
READ_LIST="References/Goodall/SraRunTable.txt" #https://www.ebi.ac.uk/ena/browser/view/PRJEB24436
READ_FOLDER="References/Goodall/Reads" #should be empty
MIN_READ_LENGTH="25"
TBL2BED="./tbl2bed.awk"
SELE="sele/sele"
SRA_QUERY="TITLE~LB"
################################################################################
# steps that only need to be set up once
################################################################################
# create conda environment, only need to do the first time
conda create \
	-n transit \
	-c bioconda \
	'bowtie=1.2.3' sra-tools entrez-direct git 'transit==3.1.0' 'samtools==1.10' 'bedtools==2.29.2' 'seqkit==0.11.0' 'exonerate==2.4.0' tabulate 'python>3'

# get a copy of my tsv tools
git clone https://github.com/ryandward/sele.git

#make sure to have csv2tsv
cpanm App::csv2tsv

################################################################################

# create read folder
mkdir $READ_FOLDER

# activate the environment
conda activate transit

# get the reads
prefetch $(csv2tsv $READ_LIST | $SELE -q -c Run -w $SRA_QUERY) --output-directory $READ_FOLDER

# dump the fastq
for i in $READ_FOLDER/*; do
	fasterq-dump --split-files $i --outdir $READ_FOLDER;
done

# get genome fasta for target
for i in $GENOME_TARGET; do
	efetch -db nuccore -format fasta -id $i > ${i}.fasta &&
	  file ${i}.fasta |
	  grep -iq ascii &&
	    echo "${i}.fasta contains data, continue to next step." ||
	    	echo "Emtpy file: ${i}.fasta, try efetch step again."
done;

# get feature table for target
for i in $GENOME_TARGET; do
	efetch -db nuccore -format ft -id $i > ${i}.tbl &&
	  file ${i}.tbl |
	  grep -iq ascii &&
	    echo "${i}.tbl contains data, continue to next step." ||
	    	echo "Emtpy file: ${i}.tbl, try efetch step again."
done;
################################################################################
for i in $GENOME_TARGET; do
	bowtie-build -o1 ${i}.fasta $i
done;

for i in $GENOME_TARGET; do
	bowtie \
		-S \
		-m2 \
		--threads 12 \
		--no-unal \
		--best \
		--strata \
		--tryhard \
		$i <(awk -v MRL=$MIN_READ_LENGTH 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= MRL) {print header, seq, qheader, qseq}}' \
		$READ_FOLDER/*.fastq) |
		samtools view -S -b -F 4 > ${i}_aln.bam
done;

# sort the bam then save it to a temp file
for i in $GENOME_TARGET; do
	samtools sort ${i}_aln.bam > ${i}_aln.bam.sort
done;

# move the temp file back to the original bam, to save space
for i in $GENOME_TARGET; do
	mv ${i}_aln.bam.sort ${i}_aln.bam
done;

################################################################################

for i in $GENOME_TARGET; do
	bam2wig.pl \
		--rpm \
	  --in ${i}_aln.bam \
	  --start \
	  --out ${i}_aln.wig
done;

################################################################################

for i in $GENOME_TARGET; do
	$TBL2BED ${i}.tbl |
	awk '$7 ~ "CDS" {print $1,$2+1,$3,$6,($3-$2)/3,$7,$8,$5,$4}' > ${i}.prot_table
done;

for i in $GENOME_TARGET; do
	$TBL2BED ${i}.tbl |
	awk '$7 ~ "CDS"' > ${i}_coding.bed
done;

for i in $GENOME_TARGET; do
	transit tn5gaps \
		-r mean \
		${i}_aln.wig \
		${i}.prot_table \
		${i}_transit_tn5gaps.tsv
done;

rm ${ORGANISM}_transit_tn5gaps.tsv

for i in $GENOME_TARGET; do
	cat	${i}_transit_tn5gaps.tsv | grep -v "#" >> ${ORGANISM}_transit_tn5gaps.tsv
done;

mkdir ${ORGANISM}_workfiles

for i in $GENOME_TARGET;
	do mv $i* ${ORGANISM}_workfiles;
done
