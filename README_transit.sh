GENOME_TARGET="NZ_CP012004.1"
MIN_READ_LENGTH="25"
READ_FOLDER="/home/ryandward/Git/transit/SRR"
TBL2BED="/home/ryandward/Git/transit/tbl2bed.awk"

conda create \
	-n transit \
	-c bioconda \
	'bowtie=1.2.3' entrez-direct git 'transit==3.1.0' 'samtools==1.10' 'bedtools==2.29.2' 'seqkit==0.11.0' 'exonerate==2.4.0' tabulate 'python>3'

conda activate transit

efetch -db nuccore -format fasta -id $GENOME_TARGET > ${GENOME_TARGET}.fasta &&
  file ${GENOME_TARGET}.fasta |
  grep -iq ascii &&
    echo "${GENOME_TARGET}.fasta contains data, continue to next step." ||
    	echo "Emtpy file: ${GENOME_TARGET}.fasta, try efetch step again."

efetch -db nuccore -format ft -id $GENOME_TARGET > ${GENOME_TARGET}.tbl &&
  file ${GENOME_TARGET}.tbl |
  grep -iq ascii &&
    echo "${GENOME_TARGET}.tbl contains data, continue to next step." ||
    	echo "Emtpy file: ${GENOME_TARGET}.tbl, try efetch step again."

$TBL2BED ${GENOME_TARGET}.tbl |
awk '$7 ~ "CDS" {print $1,$2,$3,$6,($3-$2)/3,$7,$8,$5,$4}' > ${GENOME_TARGET}.prot_table

$TBL2BED ${GENOME_TARGET}.tbl |
awk '$7 ~ "CDS"' > ${GENOME_TARGET}_coding.bed

bowtie-build -o1 ${GENOME_TARGET}.fasta $GENOME_TARGET

bowtie \
	-S \
	-m2 \
	--threads 12 \
	--no-unal \
	--best \
	--strata \
	--tryhard \
	$GENOME_TARGET <(awk -v MRL=$MIN_READ_LENGTH 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= MRL) {print header, seq, qheader, qseq}}' \
	$READ_FOLDER/*.fastq) |
	samtools view -S -b -F 4 > ${GENOME_TARGET}_aln.bam

# sort the bam then save it to a temp file
samtools sort ${GENOME_TARGET}_aln.bam > ${GENOME_TARGET}_aln.bam.sort

# move the temp file back to the original bam, to save space
mv ${GENOME_TARGET}_aln.bam.sort ${GENOME_TARGET}_aln.bam

bam2wig.pl \
  --in ${GENOME_TARGET}_aln.bam \
  --start \
  --rpm \
  --out ${GENOME_TARGET}_aln.wig

transit tn5gaps \
	-m2 \
	-r mean \
	${GENOME_TARGET}_aln.wig \
	${GENOME_TARGET}.prot_table \
	${GENOME_TARGET}_transit_tn5gaps.tsv
