#!/bin/bash
                                                                         
echo "
▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄


      ▀███▀▀▀██████▀▀██▀▀██████▀▀▀███▀███▄   ▀███▀████▀   ▀███▄   ▀███▀ ▄▄█▀▀▀█▄█
        ██    ▀██▀   ██   ▀█ ██    ▀█  ███▄    █   ██       ███▄    █ ▄██▀     ▀█
        ██   █       ██      ██   █    █ ███   █   ██       █ ███   █ ██▀       ▀
        ██████       ██      ██████    █  ▀██▄ █   ██       █  ▀██▄ █ ██         
        ██   █  ▄    ██      ██   █  ▄ █   ▀██▄█   ██     ▄ █   ▀██▄█ ██▄        
        ██     ▄█    ██      ██     ▄█ █     ███   ██    ▄█ █     ███ ▀██▄     ▄▀
      ▄██████████  ▄████▄  ▄█████████████▄    ██ █████████████▄    ██   ▀▀█████▀ 
      - Team EvolOMICS
         
     	   END TO END NOVEL LNCRNA IDENTIFICATION AND QUANTIFICATION PIPELINE	   

▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄

"
while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please specify an ID for this run. It can be a name/number or experiment/project name" 
read run

echo "You have entered $run as your run ID. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

path="/$run/"


while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter your raw reads (.fastq.gz) location" 
read fastq_query

echo "You have entered $fastq_query as your reads location. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter your lncRNAs fasta (.fa) file location" 
read lncRNAs

echo "You have entered $lncRNAs as your lncRNAs fasta location. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter your protein coding transcripts fasta (.fa) file location" 
read pc_t

echo "You have entered $pc_t as your protein coding transcripts file location. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter the location to the R script 'tximport_deseq.R'" 
read delnc

echo "You have entered $de_lnc as the location to the R script 'tximport_deseq.R'. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter the location to the experimental design file 'Samples.txt'" 
read samples

echo "You have entered $samples as the location to the design file 'Samples.txt'. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter the maximum number of threads to run this pipeline" 
read threads

echo "You have entered $threads as the maximum number of threads. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING ETENLNC"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING DIRECTORIES"

mkdir -p /user_data/$run/results/fastqc/
mkdir -p /user_data/$run/results/fastp/
mkdir -p /user_data/$run/results/fastp_reports/
mkdir -p /user_data/$run/results/DE/mRNA/
mkdir -p /user_data/$run/results/DE/lncRNA/known/
mkdir -p /user_data/$run/results/DE/lncRNA/novel/
mkdir -p /user_data/$run/results/quant/
mkdir -p /user_data/$run/results/novel_quant/
mkdir -p /user_data/$run/results/mrna_quant/
mkdir -p /user_data/$run/query/fastq/
mkdir -p /user_data/$run/query/lncRNAs/
mkdir -p /user_data/$run/query/protein_coding_transcripts/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DIRECTORIES CREATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "ASSIGNING DIRECTORIES"

fastqc=/user_data/$run/results/fastqc/
fastp=/user_data/$run/results/fastp/
fastpr=/user_data/$run/results/fastp_reports/
fastq=/user_data/$run/query/fastq/
de_mrna=/user_data/$run/results/DE/mRNA/
de_klnc=/user_data/$run/results/DE/lncRNA/known/
de_nlnc=/user_data/$run/results/DE/lncRNA/novel/
klncRNA=/user_data/$run/query/lncRNAs/
pmrna=/user_data/$run/query/protein_coding_transcripts/
base=/user_data/$run/results/quant/
base2=/user_data/$run/results/novel_quant/
base3=/user_data/$run/results/mrna_quant/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING SYMLINKS"

ln -s $fastq_query $fastq
ln -s $lncRNAs $klncRNA
ln -s $pc_t $pmrna
touch /user_data/$run/summary.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING QUALITY CHECK & PROCESSING"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING FASTP"

cd $fastq

find $fastq -name "*.fastq.gz" | sort | paste - - | while read A B

do
a=`basename ${A} | sed 's/.sra_1/_1/' | awk -F "." '{print $1}'`
b=`basename ${B} | sed 's/.sra_2/_2/' | awk -F "." '{print $1}'`

echo [`date +"%Y-%m-%d %H:%M:%S"`] "TRIMMING FILES"

fastp --thread=$threads --length_required=10 --qualified_quality_phred=32 --in1=${A} --in2=${B} --out1=$a\_trimmed.fastq.gz --out2=$b\_trimmed.fastq.gz --json=$a.json --html=$a.html

#--thread= number of worker threads (max 16)
#USE your required adapter after -a, default = automatic detection

mv -v $a\_trimmed.fastq.gz $b\_trimmed.fastq.gz $fastp

echo ""
echo "--------------------------------$a AND $b TRIMMED--------------------------------"

done
mv *.html "$fastpr"  
mv *.json "$fastpr" 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "REPORTS GENERATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING FASTQC"

cd $fastp
fastqc *.gz -t $threads
mv *.html "$fastqc"  
mv *.zip "$fastqc" 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QC REPORTS GENERATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUALITY CHECK AND PROCESSING DONE"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING QUANTIFICATION OF LNCRNAS"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING SALMON"

cd $base

echo [`date +"%Y-%m-%d %H:%M:%S"`] "INDEXING KNOWN LNCRNAS"

salmon index -t $lncRNAs -i index -k 31

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFYING...."

find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

salmon quant -i index -l A -1 ${A} -2 ${B} --validateMappings -o $a

done

#==================================================================================================================================================

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING QUANTIFICATION OF mRNAs"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING SALMON"

cd $base3

echo [`date +"%Y-%m-%d %H:%M:%S"`] "INDEXING mRNAs"

salmon index -t $pc_t -i index -k 31

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFYING...."

find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

salmon quant -i index -l A -1 ${A} -2 ${B} --validateMappings -o $a

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFICATION DONE!!!!"

#==================================================================================================================================================

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING DIFFERENTIAL EXPRESSION ANALYSIS"

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Starting differential analysis of known lncRNAs
cd $base
cp $delnc $samples $base
Rscript tximport_deseq.R $base $run vst

#Fixes left shift
mv -v "$run"\ .csv bak_"$run"\ .csv
echo -n "Gene_ID", > $run.csv
cat bak_"$run"\ .csv >> $run.csv

#Filters significant DEGs
cat $run.csv | awk -F ',' '$7<0.05{print$0}' > "$run"_significant.csv

#FIlters significant upregulated and downregulated DEGs
awk -F ',' '$3>0{print$0}' "$run"_significant.csv > "$run"_significant_up.csv
awk -F ',' '$3<0{print$0}' "$run"_significant.csv > "$run"_significant_down.csv
cat "$run"_significant_down.csv "$run"_significant_up.csv > "$run"_significant_DE.csv
cat "$run"_significant_DE.csv | awk -F ',' '{print$1}' | sed -e 's/"//g' > DE_lnc_IDs.txt
seqtk subseq $lncRNAs DE_lnc_IDs.txt > DE_lnc.fa

cp "$run"_significant_up.csv "$run"_significant_down.csv "$run"_significant.csv DE_lnc.fa $de_klnc

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Starting differential analysis of mRNAs
cp $delnc $samples $base3
cd $base3
Rscript tximport_deseq.R $base3 $run varianceStabilizingTransformation

#Fixes left shift
mv -v "$run"\ .csv bak_"$run"\ .csv
echo -n "Gene_ID", > $run.csv
cat bak_"$run"\ .csv >> $run.csv

#Filters significant DEGs
cat $run.csv | awk -F ',' '$7<0.05{print$0}' > "$run"_significant.csv

#FIlters significant upregulated and downregulated DEGs
awk -F ',' '$3>0{print$0}' "$run"_significant.csv > "$run"_significant_up.csv
awk -F ',' '$3<0{print$0}' "$run"_significant.csv > "$run"_significant_down.csv
cat "$run"_significant_down.csv "$run"_significant_up.csv > "$run"_significant_DE.csv
cat "$run"_significant_DE.csv | awk -F ',' '{print$1}' | sed -e 's/"//g' > DE_m_IDs.txt
seqtk subseq $pc_t DE_m_IDs.txt > DE_mrna.fa

cp "$run"_significant_up.csv "$run"_significant_down.csv "$run"_significant.csv DE_mrna.fa $de_mrna

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DE COMPLETE"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "FINISHING UP"

cd $base
p=`cat DE_lnc.fa | grep '>' | wc -l`
cd $base3
n=`cat DE_mrna.fa | grep '>' | wc -l`

echo "ETENLNC RUN SUMMARY" >> /user_data/$run/summary.txt
echo "PLEASE FIND BELOW THE DETAILS FOR THE PRESENT RUN" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "Current Run ID: $run" >> /user_data/$run/summary.txt
echo "Total known lncRNAs found to be differentially expressed: $p" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "Total mRNAs found to be differentially expressed: $n" >> /user_data/$run/summary.txt
echo "End_of_run" >> /user_data/$run/summary.txt

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "CLEANING UP"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "END OF PIPELINE"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"

##################################################################################################################################################
###################################################################EMD OF SCRIPT##################################################################
##################################################################################################################################################

