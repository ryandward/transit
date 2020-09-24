GENOME_TARGET_1="U00096.3"
GENOME_TARGET_2="NZ_CP012004.1"
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

##############################################################################

# convert the feature tables into bed files
$TBL2BED ${GENOME_TARGET_1}.tbl |
awk '$7 ~ "CDS"' > ${GENOME_TARGET_1}_coding.bed

$TBL2BED ${GENOME_TARGET_2}.tbl |
awk '$7 ~ "CDS"' > ${GENOME_TARGET_2}_coding.bed

# convert the bedfiles into fasta nucleotide sequences
bedtools getfasta -nameOnly -s -bed ${GENOME_TARGET_1}_coding.bed -fi ${GENOME_TARGET_1}.fasta |
sed 's/(.)//g'> ${GENOME_TARGET_1}_coding.fna

bedtools getfasta -nameOnly -s -bed ${GENOME_TARGET_2}_coding.bed -fi ${GENOME_TARGET_2}.fasta |
sed 's/(.)//g' > ${GENOME_TARGET_2}_coding.fna

# some nucleotide sequences are discontiguous, so merge them
awk 'FNR%2==1 && !($0 in seen){print; seen[$0]; getline; print} FNR%2==1 && $0 in seen {getline; print}' ${GENOME_TARGET_1}_coding.fna |
awk 'BEGIN{ORS=""} /^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' > tmp && mv tmp ${GENOME_TARGET_1}_coding.fna

awk 'FNR%2==1 && !($0 in seen){print; seen[$0]; getline; print} FNR%2==1 && $0 in seen {getline; print}' ${GENOME_TARGET_2}_coding.fna |
awk 'BEGIN{ORS=""} /^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' > tmp && mv tmp ${GENOME_TARGET_2}_coding.fna

# translate them into amino acid sequences
seqkit translate -f1 -x -w0 -T11 -x -M ${GENOME_TARGET_1}_coding.fna > ${GENOME_TARGET_1}_coding.faa

seqkit translate -f1 -x -w0 -T11 -x -M ${GENOME_TARGET_2}_coding.fna > ${GENOME_TARGET_2}_coding.faa

# use exonerate to generate a list of pair-wise best matches
# output the query, target, percent similarity, raw score
exonerate --bestn 1 ${GENOME_TARGET_1}_coding.faa ${GENOME_TARGET_2}_coding.faa --showalignment 0 --showvulgar 0 --ryo "%qi\t%ti\t%ps\t%s\t"\#ortho"\n" > ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv

exonerate --bestn 1 ${GENOME_TARGET_2}_coding.faa ${GENOME_TARGET_1}_coding.faa --showalignment 0 --showvulgar 0 --ryo "%ti\t%qi\t%ps\t%s\t"\#ortho"\n" >> ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv

# find orthos that occur more than once in the list, i.e. pair-wise best
awk '$5=="#ortho" {seen[$1 OFS $2 OFS $3 OFS $4]++} END{for(i in seen){if(seen[i]>1){print i}}}' \
  ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos.tsv \
  > ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best.tsv

awk 'NR == FNR {def[$3]=$1 OFS $2} NR != FNR && $1 in def{print $0, def[$1]}' \
  ${GENOME_TARGET_1}_genes.tsv ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best.tsv |
   sort -k6,6V \
   > ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv




############ENRICHMENT ANALYSIS#################################################

mkdir Enrichment 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos 2>/dev/null

# add gene names
awk 'NR == FNR {def[$3]=$1 OFS $2} NR != FNR && $1 in def{print $0, def[$1]}' \
  ${GENOME_TARGET_1}_genes.tsv ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best.tsv |
  sort -k6,6V \
  > ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv

# only use gene namess
awk '{print $6}' \
  ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv \
  > Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_sym.tsv

# number of genes
BACKGROUND="$(wc -l ${GENOME_TARGET_1}_coding.bed | cut -d' ' -f1)"

# find enrichment of orthologs
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_sym.tsv \
  -t=sym \
  -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos/ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

# compare essentiality of both organisms wrt the first
awk 'NR==FNR {ess[$1]=$NF} NR != FNR {print $0, ess[$1]}' \
  ${GENOME_TARGET_1}_transit_tn5gaps.tsv ${GENOME_TARGET_1}_${GENOME_TARGET_2}_orthos_best_ext.tsv |
    awk 'NR==FNR {ess[$1]=$NF} NR != FNR {print $0, ess[$2]}' \
      ${GENOME_TARGET_2}_transit_tn5gaps.tsv - |
        awk '{print $1, $7, $2, $8, $5, $6}' > Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv

# create empty directories
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TP 2>/dev/null
mkdir Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_TN 2>/dev/null


# find false positives genome 2 wrt 1
awk '$2=="Non-essential" && $4=="Essential" {print $NF}' \
  Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv > \
    Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP.tsv

# find false negatives genome 2 wrt 1
awk '$2=="Essential" && $4=="Non-essential" {print $NF}' \
  Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_essential_compare.tsv > \
    Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FN.tsv

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
~/Git/geneSCF/geneSCF-master-v1.1-p2/geneSCF-master-source-v1.1-p2/geneSCF \
  -m=normal \
  -i=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP.tsv \
  -t=sym \
  -o=Enrichment/${GENOME_TARGET_1}_${GENOME_TARGET_2}_FP/ \
  -db=GO_all \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

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
