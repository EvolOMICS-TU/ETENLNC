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

script /user_data/log.txt

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
                                                                                                                                                                                                                                                                                
echo "Please enter the strandedness of your data. Use 'RF' for reverse-stranded, 'FR' for forward-stranded and 'UN' for unstranded" 
read strand

echo "You have entered $strand as the strandedness of your data. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter your reference genome (.fa) location" 
read ref_fa

echo "You have entered $ref_fa as your reference genome location. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter your reference annotation (.gtf) location" 
read ref_gtf

echo "You have entered $ref_gtf as your reference annotation location. Do you wish to continue (Y/N)?" 
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
                                                                                                                                                                                                                                                                                
echo "Please enter your known miRNA fasta (.fa) location" 
read mirna

echo "You have entered $mirna as your known miRNA fasta (.fa) location. Do you wish to continue (Y/N)?" 
read y

case $y in 
	Y ) echo "";
	  break;;
	* ) echo "";;

esac

done

while true;
do     
                                                                                                                                                                                                                                                                                
echo "Please enter your known proteins fasta (.fa) location" 
read proteins

echo "You have entered $proteins as your known protein fasta (.fa) location. Do you wish to continue (Y/N)?" 
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

mkdir -p /user_data/$run/results/fastqc/ /user_data/$run/results/
mkdir -p /user_data/$run/results/fastp/
mkdir -p /user_data/$run/results/fastp_reports/
mkdir -p /user_data/$run/results/DE/mRNA/
mkdir -p /user_data/$run/results/DE/lncRNA/known/
mkdir -p /user_data/$run/results/DE/lncRNA/novel/
mkdir -p /user_data/$run/results/stringtie/
mkdir -p /user_data/$run/results/hisat2/
mkdir -p /user_data/$run/results/samtools/
mkdir -p /user_data/$run/results/gffcompare/
mkdir -p /user_data/$run/results/lnc_classes/
mkdir -p /user_data/$run/results/BLAST/
mkdir -p /user_data/$run/results/Predicted_LncRNAs/
mkdir -p /user_data/$run/results/index/
mkdir -p /user_data/$run/results/lists/
mkdir -p /user_data/$run/results/filters/
mkdir -p /user_data/$run/results/quant/
mkdir -p /user_data/$run/results/novel_quant/
mkdir -p /user_data/$run/results/mrna_quant/
mkdir -p /user_data/$run/results/downstream/LncTar
mkdir -p /user_data/$run/results/downstream/SEEKR/
mkdir -p /user_data/$run/results/downstream/LPI/
mkdir -p /user_data/$run/results/downstream/RNAFold/
mkdir -p /user_data/$run/results/downstream/miranda/
mkdir -p /user_data/$run/query/fastq/
mkdir -p /user_data/$run/query/reference_genome/
mkdir -p /user_data/$run/query/reference_annotation/
mkdir -p /user_data/$run/query/lncRNAs/
mkdir -p /user_data/$run/query/protein_coding_transcripts/
mkdir -p /user_data/$run/query/mirna/
mkdir -p /user_data/$run/query/proteins/
touch /user_data/$run/results/lists/list1.txt
touch /user_data/$run/results/lists/list2.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DIRECTORIES CREATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "ASSIGNING DIRECTORIES"

fastqc=/user_data/$run/results/fastqc/
fastp=/user_data/$run/results/fastp/
fastpr=/user_data/$run/results/fastp_reports/
hisat2=/user_data/$run/results/hisat2/
samtools=/user_data/$run/results/samtools/
fastq=/user_data/$run/query/fastq/
reffa=/user_data/$run/query/reference_genome/
refgtf=/user_data/$run/query/reference_annotation/
index=/user_data/$run/results/index/
de_mrna=/user_data/$run/results/DE/mRNA/
de_klnc=/user_data/$run/results/DE/lncRNA/known/
de_nlnc=/user_data/$run/results/DE/lncRNA/novel/
stringtie=/user_data/$run/results/stringtie/
list1=/user_data/$run/results/lists/list1.txt
list2=/user_data/$run/results/lists/list2.txt
gffcompare=/user_data/$run/results/gffcompare/
gff=/user_data/$run/results/gffcompare/gffcompare-0.11.4.Linux_x86_64/
classes=/user_data/$run/results/lnc_classes/
filters=/user_data/$run/results/filters/
klncRNA=/user_data/$run/query/lncRNAs/
kmirna=/user_data/$run/query/mirna/
kproteins=/user_data/$run/query/proteins/
pmrna=/user_data/$run/query/protein_coding_transcripts/
blast=/user_data/$run/results/BLAST/
novel=/user_data/$run/results/Predicted_LncRNAs/
base=/user_data/$run/results/quant/
base2=/user_data/$run/results/novel_quant/
base3=/user_data/$run/results/mrna_quant/
down=/user_data/$run/results/downstream/
SEEKR=/user_data/$run/results/downstream/SEEKR
LPI=/user_data/$run/results/downstream/LPI
RNAfold=/user_data/$run/results/downstream/RNAFold
miranda=/user_data/$run/results/downstream/miranda
LncTar=/user_data/$run/results/downstream/LncTar
filter=/script/200ntfilter.pl
gffcom=/gffcompare/gffcompare
cpc=/CPC2/bin/CPC2.py
lnctarr=/LncTar/LncTar.pl
seekrr=/seekr/seekr/kmer_counts.py
rnafold=/RNAfold
miraanda=/miranda
capsule=/LPI/main.py
lpic=/fasta_data/
delnc=/script/tximport_deseq.R

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING SYMLINKS"

ln -s $fastq_query $fastq
ln -s $ref_fa $reffa
ln -s $ref_gtf $refgtf
ln -s $lncRNAs $klncRNA
ln -s $mirna $kmirna
ln -s $proteins $kproteins
ln -s $pc_t $pmrna
touch /user_data/$run/summary.txt

cd $kproteins
mv *.fa proteins.fa

ln -s $kproteins/proteins.fa $lpic

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
echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING ALIGNMENT AND ASSEMBLY"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING HISAT2"

# =================================================================

echo [`date +"%Y-%m-%d %H:%M:%S"`] "BUILDING INDICES"

cd $index
hisat2-build -p $threads $reffa/*.fa index


# =================================================================

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING ALIGNMENT"


find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

case $strand in 
	RF ) 
	hisat2 --rna-strandness R --threads $threads --dta -x $index/index -1 ${A} -2 ${B} -S  $hisat2/$a.sam;;
	FR )
	hisat2 --rna-strandness F --threads $threads --dta -x $index/index -1 ${A} -2 ${B} -S  $hisat2/$a.sam;;
	UN ) 
	hisat2 --threads $threads --dta -x $index/index -1 ${A} -2 ${B} -S  $hisat2/$a.sam;;
	* ) 
	echo "Please check the strandedness"
	echo "HISAT2: Please check the strandedness" >> ~/$run/summary.txt;;
esac

done

#END

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CONVERTING TO BAM"

cd $hisat2

for i in *.sam
do 
echo "Converting $i"
samtools sort -@ "$threads" -o $samtools/$i.bam $i		
echo "$i converted"
done

cd $samtools

# =================================================================

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING STRINGTIE"

#cd $stringtie
for i in `ls *.bam`
do

case $strand in 
	RF ) 
	stringtie --rf -p $threads -G $ref_gtf $i -v -o $i.gtf;;
	FR ) 
	stringtie --fr -p $threads -G $ref_gtf $i -v -o $i.gtf;;
	UN ) 
	stringtie -p $threads -G $ref_gtf $i -v -o $i.gtf;;
	* ) 
	echo "Please check the strandedness"
	echo "Stringtie: Please check the strandedness" >> ~/$run/summary.txt;;
esac

echo "$i.gtf">> $list1
done

stringtie --merge -G $ref_gtf -p $threads -v -o $stringtie/stringtie_merged.gtf $list1  	#--rf-> strandedness, reverse strand; --fr-> forward strand
echo "stringtie_merged.gtf" >> $list2

cp *.bam.gtf $stringtie 
cd $stringtie
cp stringtie_merged.gtf $gffcompare

echo [`date +"%Y-%m-%d %H:%M:%S"`] "ALIGNMENT AND ASSEMBLY DONE"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING LNC IDENTIFICATION"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING GFFCOMPARE"

cd $gffcompare
$gffcom -r $ref_gtf -o $gffcompare/gffannotated.gtf -i $list2
cp gffannotated.gtf.annotated.gtf gffannotated.gtf
cp gffannotated.gtf $classes
echo [`date +"%Y-%m-%d %H:%M:%S"`] "ISOLATING CLASSES"

cd $classes 
cat gffannotated.gtf | grep 'class_code "[ioux]"' > selected.gtf
cat gffannotated.gtf | grep 'class_code "x"' > x.gtf
cat gffannotated.gtf | grep 'class_code "i"' > i.gtf
cat gffannotated.gtf | grep 'class_code "u"' > u.gtf
cat gffannotated.gtf | grep 'class_code "o"' > o.gtf

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CONVERTING TO FASTA"

cp $ref_fa $classes
mv *.fa ref.fa
#bedtools getfasta -fi ref.fa -fo selected.fa -bed selected.gtf
gffread -w selected.fa -g ref.fa selected.gtf

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FILTERING WITH CRITERIA"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FILTERING TRANSCRIPTS GREATER THAN 200 NUCLEOTIDES"

cd $filters
cp $filter $filters
cp $classes/selected.fa $filters
cp $gffcompare/gffannotated.gtf $filters
cp $classes/ref.fa $filters
perl $filter 200 selected.fa > gt200.fa

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING CODING POTENTIAL ANALYSIS"

cd $filters

python3 $cpc -i gt200.fa -o CPC
cat CPC.txt | grep 'noncoding' > noncodingRNA_1.txt
cat noncodingRNA_1.txt | awk '$3 < 300 {print$0}' > noncodingRNA.txt
cat noncodingRNA.txt | awk '{print $1}' > noncodeIDs.txt
sed 's/^/transcript_id "&/g' noncodeIDs.txt > noncode_2.txt
cat gffannotated.gtf | fgrep -f noncode_2.txt > Noncoding.gtf

cat Noncoding.gtf | awk '$13 == "exon_number" {print $14, $10}' > exon1.txt
cat exon1.txt | awk -F',' '{gsub(/"/, "", $1); print $0}' | awk '$1 > 2 {print $0}' > exon2.txt
cat exon2.txt | awk -F ';' '{print $2}' | sed 's/ //g' > exon3.txt
sed 's/^/transcript_id "&/g' exon3.txt > exon4.txt
cat gffannotated.gtf | fgrep -f exon4.txt > Exon_filtered.gtf
gffread -w Noncoding200nt.fa -g ref.fa Exon_filtered.gtf

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FILTERING DONE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING BLAST AGAINST KNOWN LNCRNAS"

cd $blast
cp $klncRNA/*.fa $blast
mv *.fa lncRNAs.fasta
cp $filters/Noncoding200nt.fa $blast
cp $filters/exon3.txt $blast
cp $filters/ref.fa $blast
cp $filters/gffannotated.gtf $blast

makeblastdb -in Noncoding200nt.fa -input_type fasta -parse_seqids -dbtype nucl -out LncRNA_novel
blastn -db LncRNA_novel -query $lncRNAs -out BLAST.txt -evalue 0.001 -outfmt 6 -word_size 7 -num_threads $threads
cat BLAST.txt | awk '{print $2}' > BLASTids.txt
cat BLASTids.txt | awk '!seen[$0]++' > BLAST_nr.txt
grep -Fvx -f BLAST_nr.txt exon3.txt > novel_lncRNAs.txt
sed 's/^/transcript_id "&/g' novel_lncRNAs.txt > LncRNA_2.txt
cat gffannotated.gtf | fgrep -f LncRNA_2.txt > Novel_LncRNAs.gtf
gffread -w Novel_LncRNAs.fa -g ref.fa Novel_LncRNAs.gtf
#bedtools getfasta -fi ref.fa -fo Novel_LncRNAs.fa -bed Novel_LncRNAs.gtf

sed 's/^/transcript_id "&/g' BLAST_nr.txt > BLAST_inter1.txt
cat gffannotated.gtf | fgrep -f BLAST_inter1.txt > BLAST_nr.gtf
gffread -w Known_LncRNAs.fa -g ref.fa BLAST_nr.gtf
#bedtools getfasta -fi ref.fa -fo Known_LncRNAs.fa -bed BLAST_nr.gtf

cp Novel_LncRNAs.fa $novel
cp Novel_LncRNAs.fa $lpic

echo [`date +"%Y-%m-%d %H:%M:%S"`] "PUTATIVE LNCRNAS IDENTIFIED"

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

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING QUANTIFICATION OF NOVEL LNCRNAS"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING SALMON"

cd $base2

echo [`date +"%Y-%m-%d %H:%M:%S"`] "INDEXING NOVEL LNCRNAs"

salmon index -t $novel/Novel_LncRNAs.fa -i index -k 31

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFYING...."

find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

salmon quant -i index -l A -1 ${A} -2 ${B} --validateMappings -o $a

done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFICATION DONE!!!!"

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
cat "$run"_significant_DE.csv | awk -F ',' '{print$1}' | sed -e 's/"//g' > DE_klnc_IDs.txt
seqtk subseq $lncRNAs DE_klnc_IDs.txt > DE_klnc.fa

cp "$run"_significant_up.csv "$run"_significant_down.csv "$run"_significant.csv DE_klnc.fa $de_klnc

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Starting differential analysis of novel lncRNAs
cp $delnc $samples $base2
cd $base2
Rscript tximport_deseq.R $base2 $run varianceStabilizingTransformation

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
cat "$run"_significant_DE.csv | awk -F ',' '{print$1}' | sed -e 's/"//g' > DE_nlnc_IDs.txt
seqtk subseq $novel/Novel_LncRNAs.fa DE_nlnc_IDs.txt > DE_nlnc.fa

cp "$run"_significant_up.csv "$run"_significant_down.csv "$run"_significant.csv DE_nlnc.fa $de_nlnc

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
seqtk subseq $pc_t DE_m_IDs.txt > DE_mrna1.fa
cut -d ' ' -f 1 DE_mrna1.fa > DE_mrna.fa

cp "$run"_significant_up.csv "$run"_significant_down.csv "$run"_significant.csv DE_mrna.fa $de_mrna

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DE COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING DOWNSTREAM ANALYSES"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DOWNSTREAM ANALYSES"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING FUNCTIONAL ANNOTATION BY SEEKR"
cd $SEEKR
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_1.csv -k 1
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_2.csv -k 2
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_3.csv -k 3
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_4.csv -k 4
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_5.csv -k 5
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_6.csv -k 6
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_7.csv -k 7
seekr_kmer_counts $de_nlnc/DE_nlnc.fa -o kmer_counts_8.csv -k 8
echo [`date +"%Y-%m-%d %H:%M:%S"`] "FUNCTIONAL ANNOTATION BY SEEKR COMPLETE"

cd $LPI
echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING PROTEIN INTERACTION PREDICTION BY CAPSULE-LPI"
cd /LPI/
python3 main.py
cp $lpic/lnc_protein.csv $LPI

#Filtering LPI results and creating nodelists
cat lnc_protein.csv | awk -F, '$4==1.0 {print $2,$3}' > lnc_protein_nodes.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "PROTEIN INTERACTION PREDICTION BY CAPSULE-LPI COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING RNA FOLD PREDICTION BY RNAFOLD"
cd $RNAfold
$rnafold < $de_nlnc/DE_nlnc.fa
echo [`date +"%Y-%m-%d %H:%M:%S"`] "RNA FOLD PREDICTION BY RNAFOLD COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING LNCRNA-miRNA PREDICTION BY miRanda"
cd $miranda
$miraanda $mirna $de_nlnc/DE_nlnc.fa -out DEnlncRNA_miRNA.txt
$miraanda $mirna $de_klnc/DE_klnc.fa -out DEklncRNA_miRNA.txt
$miraanda $mirna $de_mrna/DE_mrna.fa -out DEmRNA_miRNA.txt

#Filtering miranda results and creating nodelists
cat DEnlncRNA_miRNA.txt | grep -A 1 "Scores for this hit:" | sort | grep '>' | awk '($3>150) && ($4 < -7) {print$0}' | awk '{print $1,$2}' | cut -c 2- > mi_nlnc_nodes.txt
cat DEklncRNA_miRNA.txt | grep -A 1 "Scores for this hit:" | sort | grep '>' | awk '($3>150) && ($4 < -7) {print$0}' | awk '{print $1,$2}' | cut -c 2- > mi_klnc_nodes.txt
cat DEmRNA_miRNA.txt | grep -A 1 "Scores for this hit:" | sort | grep '>' | awk '($3>150) && ($4 < -7) {print$0}' | awk '{print $1,$2}' | cut -c 2- > mi_mrna_nodes.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "LNCRNA-miRNA PREDICTION BY miRanda COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING LNCRNA TARGET PREDICTION BY LncTar"
cd $LncTar
$lnctarr -p 1 -l $de_nlnc/DE_nlnc.fa -m $de_mrna/DE_mrna.fa -d -0.1 -s F -o Novel_targets.txt
$lnctarr -p 1 -l $de_klnc/DE_klnc.fa -m $de_mrna/DE_mrna.fa -d -0.1 -s F -o Known_targets.txt

#Filtering LncTar results and creating nodelists
cat Novel_targets.txt | awk '{print $1,$3}' > Nlnc_mRNA_nodes.txt
cat Known_targets.txt | awk '{print $1,$3}' > klnc_mRNA_nodes.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "TARGET PREDICTION COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DOWNSTREAM ANALYSES COMPLETE"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FINISHING UP"

cd $novel
p=`cat Novel_LncRNAs.fa | grep '>' | wc -l`
cd $blast
n=`cat $filters/Noncoding200nt.fa | grep '>' | wc -l`
m=`cat $lncRNAs | grep '>' | wc -l`
cd $base
a=`cat DE_klnc.fa | grep '>' | wc -l`
cd $base2
b=`cat DE_nlnc.fa | grep '>' | wc -l`
cd $base3
c=`cat DE_mrna.fa | grep '>' | wc -l`
echo "Putative LncRNAs Identified: $p"
echo "LncRNAs captured from known LncRNAs:" $(($n - $p)) "out of" $m

echo "ETENLNC RUN SUMMARY" >> /user_data/$run/summary.txt
echo "PLEASE FIND BELOW THE DETAILS FOR THE PRESENT RUN" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "Current Run ID: $run" >> /user_data/$run/summary.txt
echo "NOVEL LncRNAs" >> /user_data/$run/summary.txt
echo "Total novel lncRNAs identified in the current run: $p" >> /user_data/$run/summary.txt
echo "Total novel lncRNAs found to be differentially expressed: $b" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "KNOWN LncRNAs" >> /user_data/$run/summary.txt
echo "Total known lncRNAs captured in the current run:" $(($n - $p)) "out of" $m >> /user_data/$run/summary.txt
echo "Total known lncRNAs found to be differentially expressed: $a" >> /user_data/$run/summary.txt
echo "mRNAs" >> /user_data/$run/summary.txt
echo "Total mRNAs found to be differentially expressed: $c" >> /user_data/$run/summary.txt
echo "" >> /user_data/$run/summary.txt
echo "End_of_run" >> /user_data/$run/summary.txt

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "CLEANING UP"

rm -rf $lpic/Novel_LncRNAs.fa
rm -rf $lpic/proteins.fa

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "END OF PIPELINE"

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"
exit
##################################################################################################################################################
###################################################################EMD OF SCRIPT##################################################################
##################################################################################################################################################

