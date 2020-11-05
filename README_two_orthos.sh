ORGANISM_1="Escherichia_coli_BW25113"
ORGANISM_2="Acinetobacter_baumannii_ATCC_17978"

GENOME_TARGET_1=(CP009273.1) #https://www.ncbi.nlm.nih.gov/nuccore/CP009273.1
GENOME_TARGET_2=(CP053098.1 CP053099.1 CP053100.1) #https://www.ncbi.nlm.nih.gov/assembly/GCF_013372085.1

TBL2BED="/home/ryandward/Git/transit/tbl2bed.awk"

esearch \
  -db gene \
  -query "$GENOME_TARGET_1\[NUCL_ACCN\]" |
  efetch -format tabular |
  awk 'FNR==1{print $3, $6, $7; seen[$3]} !($3 in seen) {seen[$3]; split($7, syn, ", "); for(i in syn) {print $3, $6, syn[i], $5;}}' \
  > ${GENOME_TARGET_1}_genes.tsv

sleep 10

# make sure the one that is more likely to fail is the second one
esearch \
  -db gene \
  -query "$GENOME_TARGET_2\[NUCL_ACCN\]" |
  efetch -format tabular |
  awk 'FNR==1{print $3, $6, $7; seen[$3]} !($3 in seen) {seen[$3]; split($7, syn, ", "); for(i in syn) {print $3, $6, syn[i], $5;}}' \
  > ${GENOME_TARGET_2}_genes.tsv

for i in $GENOME_TARGET_1; do
  efetch -db nuccore -format fasta -id $i > ${i}.fasta;
done

for i in $GENOME_TARGET_2; do
  efetch -db nuccore -format fasta -id $i > ${i}.fasta;
done

##############################################################################

# convert the feature tables into bed files
# $TBL2BED ${GENOME_TARGET_1}.tbl |
# awk '$7 ~ "CDS"' > ${GENOME_TARGET_1}_coding.bed
#
# $TBL2BED ${GENOME_TARGET_2}.tbl |
# awk '$7 ~ "CDS"' > ${GENOME_TARGET_2}_coding.bed

# convert the bedfiles into fasta nucleotide sequences
for i in $GENOME_TARGET_1; do
  bedtools getfasta -nameOnly -s -bed ${i}_coding.bed -fi ${i}.fasta |
  sed 's/(.)//g'> ${i}_coding.fna
done

for i in $GENOME_TARGET_2; do
  bedtools getfasta -nameOnly -s -bed ${i}_coding.bed -fi ${i}.fasta |
  sed 's/(.)//g'> ${i}_coding.fna
done


# some nucleotide sequences are discontiguous, so merge them
for i in $GENOME_TARGET_1; do
  awk 'FNR%2==1 && !($0 in seen){print; seen[$0]; getline; print} FNR%2==1 && $0 in seen {getline; print}' ${i}_coding.fna |
  awk 'BEGIN{ORS=""} /^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' > tmp && mv tmp ${i}_coding.fna
done

for i in $GENOME_TARGET_2; do
  awk 'FNR%2==1 && !($0 in seen){print; seen[$0]; getline; print} FNR%2==1 && $0 in seen {getline; print}' ${i}_coding.fna |
  awk 'BEGIN{ORS=""} /^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' > tmp && mv tmp ${i}_coding.fna
done

# translate them into amino acid sequences
rm ${ORGANISM_1}_coding.faa
for i in $GENOME_TARGET_1; do
  seqkit translate -f1 -x -w0 -T11 -x -M ${i}_coding.fna >> ${ORGANISM_1}_coding.faa
done;

rm ${ORGANISM_2}_coding.faa
for i in $GENOME_TARGET_2; do
  seqkit translate -f1 -x -w0 -T11 -x -M ${i}_coding.fna >> ${ORGANISM_2}_coding.faa
done;

################################################################################

# use exonerate to generate a list of pair-wise best matches
# output the query, target, percent similarity, raw score
exonerate --bestn 1 ${ORGANISM_1}_coding.faa ${ORGANISM_2}_coding.faa --showalignment 0 --showvulgar 0 --ryo "%qi\t%ti\t%ps\t%s\t"\#ortho"\n" > ${ORGANISM_1}_${ORGANISM_2}_orthos.tsv

exonerate --bestn 1 ${ORGANISM_2}_coding.faa ${ORGANISM_1}_coding.faa --showalignment 0 --showvulgar 0 --ryo "%ti\t%qi\t%ps\t%s\t"\#ortho"\n" >> ${ORGANISM_1}_${ORGANISM_2}_orthos.tsv

# find orthos that occur more than once in the list, i.e. pair-wise best
awk '$5=="#ortho" {seen[$1 OFS $2 OFS $3 OFS $4]++} END{for(i in seen){if(seen[i]>1){print i}}}' \
  ${ORGANISM_1}_${ORGANISM_2}_orthos.tsv \
  > ${ORGANISM_1}_${ORGANISM_2}_orthos_best.tsv

################################################################################

rm ${ORGANISM_1}_coding.bed
for i in $GENOME_TARGET_1; do
  cat ${i}_coding.bed >> ${ORGANISM_1}_coding.bed;
done

rm ${ORGANISM_2}_coding.bed
for i in $GENOME_TARGET_2; do
  cat ${i}_coding.bed >> ${ORGANISM_2}_coding.bed;
done

################################################################################

awk 'NR == FNR {def[$4]=$5} NR != FNR && $1 in def{print $1,$2,$3,$4,def[$1]}' \
  ${ORGANISM_1}_coding.bed ${ORGANISM_1}_${ORGANISM_2}_orthos_best.tsv |
   sort -k5,5V \
   > tmp && mv tmp ${ORGANISM_1}_${ORGANISM_2}_orthos_best.tsv

############ENRICHMENT ANALYSIS#################################################

# number of genes


# compare essentiality of both organisms wrt the first
echo "Symbol\tLocus_1\tEssential_1\tP_1\tLocus_2\tEssential_2\tP_2" > ${ORGANISM_1}_${ORGANISM_2}_essentiality.tsv
awk 'NR==FNR {ess[$1]=$NF OFS $(NF-1)} NR != FNR {print $5, $1, ess[$1], $2}' \
  ${ORGANISM_1}_transit_tn5gaps.tsv ${ORGANISM_1}_${ORGANISM_2}_orthos_best.tsv |
awk 'NR==FNR {ess[$1]=$NF OFS $(NF-1)} NR != FNR {print $0, ess[$NF]}' \
  ${ORGANISM_2}_transit_tn5gaps.tsv - >> ${ORGANISM_1}_${ORGANISM_2}_essentiality.tsv

################################################################################

sele -q -c Symbol -w "Essential_1==Non-essential,Essential_2==Essential" ${ORGANISM_1}_${ORGANISM_2}_essentiality.tsv \
  > ${ORGANISM_1}_${ORGANISM_2}_diff_ess.tsv

BACKGROUND="$(sele -q ${ORGANISM_1}_coding.bed | wc -l | cut -d' ' -f1)"

~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=${ORGANISM_1}_${ORGANISM_2}_diff_ess.tsv.tsv \
  -t=sym \
  -o=./ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

################################################################################

sele -q -c Symbol -w "Essential_1==Essential,Essential_2==Non-essential" ${ORGANISM_1}_${ORGANISM_2}_essentiality.tsv \
  > ${ORGANISM_1}_${ORGANISM_2}_diff_noness.tsv

BACKGROUND="$(sele -q ${ORGANISM_2}_coding.bed | wc -l | cut -d' ' -f1)"

~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=${ORGANISM_1}_${ORGANISM_2}_diff_noness.tsv \
  -t=sym \
  -o=./ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

################################################################################

# find true positives genome 2 wrt 1
awk '$2=="Essential" && $4=="Essential" {print $NF}' \
  Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv > \
    Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TP.tsv

# find true negatives genome 2 wrt 1
awk '$2=="Non-essential" && $4=="Non-essential" {print $NF}' \
  Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv > \
    Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TN.tsv

# number of orthologs
BACKGROUND="$(wc -l ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv |
  cut -d' ' -f1)"

# find false positive enrichment


# find false negative enrichment
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN.tsv \
  -t=sym -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN/ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

# find true positive enrichment
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TP.tsv \
  -t=sym \
  -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TP/ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

# find true negative enrichment
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TN.tsv \
  -t=sym \
  -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TN/ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc


conda deactivate
