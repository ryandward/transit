

GENOME_TARGET="U00096.3"

conda create -n transit -c bioconda 'bowtie=1.2.3' entrez-direct git transit 'python>3'

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

bowtie-build -o2 ${GENOME_TARGET}.fasta $GENOME_TARGET

bowtie -S -m2 --threads 12 --no-unal --best --strata --tryhard \
$GENOME_TARGET <(cat ../Reads/ERR2249109.fastq ../Reads/ERR2249110.fastq) |
samtools view -S -b > ${GENOME_TARGET}_aln.bam

samtools sort ${GENOME_TARGET}_aln.bam > ${GENOME_TARGET}_aln.bam.sort
mv ${GENOME_TARGET}_aln.bam.sort ${GENOME_TARGET}_aln.bam

bam2wig.pl  --in ${GENOME_TARGET}_aln.bam --start --rpm --out ${GENOME_TARGET}_aln.wig

~/SamsungUSB/Scripts/tbl2bed.awk ${GENOME_TARGET}.tbl |
awk '$7=="CDS" {print $1,$2,$3,$6,($3-$2)/3,$7,$8,$5,$4}' > ${GENOME_TARGET}_aln.wig

transit tn5gaps -m2 -r mean ${GENOME_TARGET}_aln.wig ${GENOME_TARGET}.prot_table ${GENOME_TARGET}_transit_tn5gaps.tsv
