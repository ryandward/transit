cd Orthos

for i in *bed; do
  filename="${i##*/}"
  extension="${filename##*.}"
  filename="${filename%.*}"
  cat $i |
  bedtools getfasta -nameOnly -s -bed stdin -fi ${filename%.bed}.fasta \
  > "${filename%.*}".fna;
done

for i in *fna; do
  filename="${i##*/}"
  extension="${filename##*.}"
  filename="${filename%.*}"
  awk 'BEGIN{
     PROCINFO["sorted_in"] = "@ind_num_asc";
     RS="\n>";FS="\n"
   }{
     gsub (">", "", $1);
     n = split($1, parts, /[^A-z0-9_>]/);
     coords[parts[1]][NR]=$2}
     END {
       for(i in coords) {
         printf("%s%s\n", ">", i);
         for(j in coords[i]) {
           printf ("%s", coords[i][j])
         }
         printf ("\n");
       }
     }' $i \
      > "${filename%.*}".tmp && mv "${filename%.*}".tmp $i;
    done

mkdir AminoAcid
for i in *fna; do
  filename="${i##*/}"
  extension="${filename##*.}"
  filename="${filename%.*}"
  seqkit translate -x -w0 -T11 -x -M $i \
    > AminoAcid/"${filename%.*}".faa ;
done

ulimit -n 524288
conda activate orthofinder &&  orthofinder -t 12 -a 12 -f ./AminoAcid

mkdir Results

ORTHOGROUPS="/home/ryandward/Class/Genetics885/transit/Orthos/AminoAcid/Results_Dec02/Orthogroups.csv"
dos2unix $ORTHOGROUPS
UNASSIGNED="/home/ryandward/Class/Genetics885/transit/Orthos/AminoAcid/Results_Dec02/Orthogroups_UnassignedGenes.csv"
dos2unix $UNASSIGNED

awk -vFS='\t' \
    'NR == 1{
      for(i=2;i<=NF;i++) {
        species[i]=$i;
      }
      next;
    }
    NR>1{
      for(i in species) {
        split($i, genes, ", ");
        for(j in genes) {
          print $1, species[i], genes[j];
        }
      }
    }' $ORTHOGROUPS $UNASSIGNED > Results/linear_orthos.tsv

cd Results

cat linear_orthos.tsv | awk '{orthos[$2][$1]++} END{for(i in orthos){for(j in orthos[i]){print j,i,orthos[i][j]}}}' | sort -V | awk 'NR==FNR{orthos[$1][$2]=$3} NR!=FNR{print $0, orthos[$1][$2]}' - linear_orthos.tsv | awk '$1!=""' | sort -V > linear_orthos_paralogs.tsv
cat *tn5gaps.tsv | awk 'NR==FNR{name[$1]=$2; ess[$1]=$10} FNR!=NR{print $0, name[$3], ess[$3]}' - linear_orthos_paralogs.tsv > linear_orthos_paralogs_score.tsv

FDR=0.05
STRAIN="Escherichia_coli_BW25113"
BACKGROUND=$(awk -v STRAIN=$STRAIN '$2==STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk -v FDR=$FDR -v STRAIN=$STRAIN \
	'$2=="Escherichia_coli_BW25113" && $6<FDR && $4==1 {print $5}' \
	linear_orthos_paralogs_score.tsv > ${STRAIN}_FDR${FDR}

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${STRAIN}_FDR${FDR} \
  -o=./ \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

#####
FDR=0.05
REF_STRAIN="Escherichia_coli_BW25113"
STRAIN="Acinetobacter_baumannii_ATCC_17978"
BACKGROUND=$(awk -v STRAIN=$STRAIN '$2==STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk '$4==1 {count[$1]++} END{for(i in count){if(count[i]==2)print i}}' linear_orthos_paralogs_score.tsv |
awk -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{valid[$1]} NR!=FNR && $1 in valid {if($2==REF_STRAIN){symbol[$1]=$5}} END{for(i in symbol){print i, symbol[i]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR -v STRAIN=$STRAIN -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{symbol[$1]=$2} NR!=FNR && $1 in symbol && $2==STRAIN && $6<FDR {$5=symbol[$1]; print $5}' - linear_orthos_paralogs_score.tsv > ${STRAIN}_FDR${FDR}

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${STRAIN}_FDR${FDR} \
  -o=./ \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

#####
FDR=0.05
REF_STRAIN="Escherichia_coli_BW25113"
STRAIN="Acinetobacter_baumannii_ATCC_17978"
BACKGROUND=$(awk -v REF_STRAIN=$REF_STRAIN '$2==REF_STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk '$4==1 {count[$1]++} END{for(i in count){if(count[i]==2)print i}}' linear_orthos_paralogs_score.tsv |
awk -v REF_STRAIN=$REF_STRAIN \
  'NR==FNR{valid[$1]} NR!=FNR && $1 in valid {if($2==REF_STRAIN){symbol[$1]=$5}} END{for(i in symbol){print i, symbol[i]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR -v STRAIN=$STRAIN -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{symbol[$1]=$2} NR!=FNR && $1 in symbol {score[$1][$2]=$6} END {for(i in symbol){print symbol[i], score[i][REF_STRAIN], score[i][STRAIN]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR '$2<FDR && $3>=FDR {print $1}' > ${REF_STRAIN}_FDR${FDR}_DIFF

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${REF_STRAIN}_FDR${FDR}_DIFF \
  -o=./ \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

#####
FDR=0.05
REF_STRAIN="Escherichia_coli_BW25113"
STRAIN="Acinetobacter_baumannii_ATCC_17978"
BACKGROUND=$(awk -v STRAIN=$STRAIN '$2==STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk '$4==1 {count[$1]++} END{for(i in count){if(count[i]==2)print i}}' linear_orthos_paralogs_score.tsv |
awk -v REF_STRAIN=$REF_STRAIN \
  'NR==FNR{valid[$1]} NR!=FNR && $1 in valid {if($2==REF_STRAIN){symbol[$1]=$5}} END{for(i in symbol){print i, symbol[i]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR -v STRAIN=$STRAIN -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{symbol[$1]=$2} NR!=FNR && $1 in symbol {score[$1][$2]=$6} END {for(i in symbol){print symbol[i], score[i][REF_STRAIN], score[i][STRAIN]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR '$2>=FDR && $3<FDR {print $1}' > ${STRAIN}_FDR${FDR}_DIFF

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${STRAIN}_FDR${FDR}_DIFF \
  -o=./ \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=ecocyc

###
cat ../../Acinetobacter_baumannii_ATCC_17978_workfiles/CP053098.1.tbl ../../Acinetobacter_baumannii_ATCC_17978_workfiles/CP053099.1.tbl ../../Acinetobacter_baumannii_ATCC_17978_workfiles/CP053100.1.tbl ../../Escherichia_coli_BW25113_workfiles/CP009273.1.tbl | ../../tbl2bed.awk | awk '$7=="CDS"' > all.bed



###
try to find unfound processes
cat unfound.txt | while read line; do grep $line ecocyc.gaf | awk -v search=$line '$9=="P"{go[$5][$3]; go[$5][search]} END{for(i in go){for (j in go[i]){print i,j}}}'; done > go_unfound.tsv
awk '{split($1,go,"~"); split($2,genes,","); for(i in genes){if(genes[i]!=""){print go[1], genes[i]}}}' /home/ryandward/Git/GeneSCF/class/lib/db/ecocyc/GO_BP_sym.txt > go_from_genescf.tsv
awk 'NR==FNR && $4~"BW25113"{need[$5]} NR!=FNR && $2 in need' all.bed go_* | sort | uniq | awk 'go[$1]{go[$1]=go[$1]","$2} !go[$1]{go[$1]=$2} END{for(i in go){print i, go[i]}}' |  awk 'NR==FNR{go[$5]=$10;} NR!=FNR{print $1"~"go[$1], $2 }' ecocyc.gaf  - | sort -V > GO_BP_sym.txt
cp GO_BP_sym.txt /home/ryandward/Git/GeneSCF/class/lib/db/curated

###
FDR=0.05
STRAIN="Escherichia_coli_BW25113"
BACKGROUND=$(awk -v STRAIN=$STRAIN '$2==STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk -v FDR=$FDR -v STRAIN=$STRAIN \
	'$2=="Escherichia_coli_BW25113" && $6<FDR && $4==1 {print $5}' \
	linear_orthos_paralogs_score.tsv > ${STRAIN}_FDR${FDR}

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${STRAIN}_FDR${FDR} \
  -o=./curated \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=curated

#####
FDR=0.05
REF_STRAIN="Escherichia_coli_BW25113"
STRAIN="Acinetobacter_baumannii_ATCC_17978"
BACKGROUND=$(awk -v STRAIN=$STRAIN '$2==STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk '$4==1 {count[$1]++} END{for(i in count){if(count[i]==2)print i}}' linear_orthos_paralogs_score.tsv |
awk -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{valid[$1]} NR!=FNR && $1 in valid {if($2==REF_STRAIN){symbol[$1]=$5}} END{for(i in symbol){print i, symbol[i]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR -v STRAIN=$STRAIN -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{symbol[$1]=$2} NR!=FNR && $1 in symbol && $2==STRAIN && $6<FDR {$5=symbol[$1]; print $5}' - linear_orthos_paralogs_score.tsv > ${STRAIN}_FDR${FDR}

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${STRAIN}_FDR${FDR} \
  -o=./curated \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=curated

#####
FDR=0.05
REF_STRAIN="Escherichia_coli_BW25113"
STRAIN="Acinetobacter_baumannii_ATCC_17978"
BACKGROUND=$(awk -v REF_STRAIN=$REF_STRAIN '$2==REF_STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk '$4==1 {count[$1]++} END{for(i in count){if(count[i]==2)print i}}' linear_orthos_paralogs_score.tsv |
awk -v REF_STRAIN=$REF_STRAIN \
  'NR==FNR{valid[$1]} NR!=FNR && $1 in valid {if($2==REF_STRAIN){symbol[$1]=$5}} END{for(i in symbol){print i, symbol[i]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR -v STRAIN=$STRAIN -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{symbol[$1]=$2} NR!=FNR && $1 in symbol {score[$1][$2]=$6} END {for(i in symbol){print symbol[i], score[i][REF_STRAIN], score[i][STRAIN]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR '$2<FDR && $3>=FDR {print $1}' > ${REF_STRAIN}_FDR${FDR}_DIFF

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${REF_STRAIN}_FDR${FDR}_DIFF \
  -o=./curated \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=curated

#####
FDR=0.05
REF_STRAIN="Escherichia_coli_BW25113"
STRAIN="Acinetobacter_baumannii_ATCC_17978"
BACKGROUND=$(awk -v STRAIN=$STRAIN '$2==STRAIN' linear_orthos_paralogs_score.tsv | wc -l)

awk '$4==1 {count[$1]++} END{for(i in count){if(count[i]==2)print i}}' linear_orthos_paralogs_score.tsv |
awk -v REF_STRAIN=$REF_STRAIN \
  'NR==FNR{valid[$1]} NR!=FNR && $1 in valid {if($2==REF_STRAIN){symbol[$1]=$5}} END{for(i in symbol){print i, symbol[i]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR -v STRAIN=$STRAIN -v REF_STRAIN=$REF_STRAIN \
	'NR==FNR{symbol[$1]=$2} NR!=FNR && $1 in symbol {score[$1][$2]=$6} END {for(i in symbol){print symbol[i], score[i][REF_STRAIN], score[i][STRAIN]}}' - linear_orthos_paralogs_score.tsv |
awk -v FDR=$FDR '$2>=FDR && $3<FDR {print $1}' > ${STRAIN}_FDR${FDR}_DIFF

~/Git/GeneSCF/geneSCF -m=normal \
  -i=${STRAIN}_FDR${FDR}_DIFF \
  -o=./curated \
  -t=sym \
  -db=GO_BP \
  -p=yes \
  -bg=$BACKGROUND \
  -org=curated
