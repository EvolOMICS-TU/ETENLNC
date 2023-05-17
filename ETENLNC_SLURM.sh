#!/bin/bash

#SBATCH --job-name=JOB1
#SBATCH --output=JOB1.out
#SBATCH --ntasks=50
                                                                    
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

#Please populate the following variables with their respective local paths

#Please specify an ID for this run. It can be a name/number or experiment/project name
run=
path="/$run/"

script /$run/log.txt

#Please enter your raw reads (.fastq.gz) location
fastq_query=

#Please enter the strandedness of your data. Use 'RF' for reverse-stranded, 'FR' for forward-stranded and 'UN' for unstranded
strand=

#Please enter your reference genome (.fa) location
ref_fa=

#Please enter your reference annotation (.gtf) location
ref_gtf=

#Please enter your lncRNAs fasta (.fa) file location
lncRNAs=

#Please enter your protein coding transcripts fasta (.fa) file location
pc_t=

#Please enter the location to the experimental design file 'Samples.txt'
samples=

#Please enter your known miRNA fasta (.fa) location
mirna=

#Please enter your known proteins fasta (.fa) location
proteins=

#Please enter the maximum number of threads to run this pipeline
threads= 

echo "▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄
"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING ETENLNC"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING DIRECTORIES"

mkdir -p ~/$run/results/fastqc/ ~/$run/results/
mkdir -p ~/$run/results/fastp/
mkdir -p ~/$run/results/fastp_reports/
mkdir -p ~/$run/results/DE/mRNA/
mkdir -p ~/$run/results/DE/lncRNA/known/
mkdir -p ~/$run/results/DE/lncRNA/novel/
mkdir -p ~/$run/results/stringtie/
mkdir -p ~/$run/results/hisat2/
mkdir -p ~/$run/results/samtools/
mkdir -p ~/$run/results/gffcompare/
mkdir -p ~/$run/results/lnc_classes/
mkdir -p ~/$run/results/BLAST/
mkdir -p ~/$run/results/Predicted_LncRNAs/
mkdir -p ~/$run/results/index/
mkdir -p ~/$run/results/lists/
mkdir -p ~/$run/results/filters/
mkdir -p ~/$run/results/quant/
mkdir -p ~/$run/results/novel_quant/
mkdir -p ~/$run/results/mrna_quant/
mkdir -p ~/$run/results/downstream/LncTar
mkdir -p ~/$run/results/downstream/SEEKR/
mkdir -p ~/$run/results/downstream/LPI/
mkdir -p ~/$run/results/downstream/RNAFold/
mkdir -p ~/$run/results/downstream/miranda/
mkdir -p ~/$run/query/fastq/
mkdir -p ~/$run/query/reference_genome/
mkdir -p ~/$run/query/reference_annotation/
mkdir -p ~/$run/query/lncRNAs/
mkdir -p ~/$run/query/protein_coding_transcripts/
mkdir -p ~/$run/query/mirna/
mkdir -p ~/$run/query/proteins/
touch ~/$run/results/lists/list1.txt
touch ~/$run/results/lists/list2.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DIRECTORIES CREATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "ASSIGNING DIRECTORIES"

fastqc=~/$run/results/fastqc/
fastp=~/$run/results/fastp/
fastpr=~/$run/results/fastp_reports/
hisat2=~/$run/results/hisat2/
samtools=~/$run/results/samtools/
fastq=~/$run/query/fastq/
reffa=~/$run/query/reference_genome/
refgtf=~/$run/query/reference_annotation/
index=~/$run/results/index/
de_mrna=~/$run/results/DE/mRNA/
de_klnc=~/$run/results/DE/lncRNA/known/
de_nlnc=~/$run/results/DE/lncRNA/novel/
stringtie=~/$run/results/stringtie/
list1=~/$run/results/lists/list1.txt
list2=~/$run/results/lists/list2.txt
gffcompare=~/$run/results/gffcompare/
gff=~/$run/results/gffcompare/gffcompare-0.11.4.Linux_x86_64/
classes=~/$run/results/lnc_classes/
filters=~/$run/results/filters/
klncRNA=~/$run/query/lncRNAs/
kmirna=~/$run/query/mirna/
kproteins=~/$run/query/proteins/
pmrna=~/$run/query/protein_coding_transcripts/
blast=~/$run/results/BLAST/
novel=~/$run/results/Predicted_LncRNAs/
base=~/$run/results/quant/
base2=~/$run/results/novel_quant/
base3=~/$run/results/mrna_quant/
down=~/$run/results/downstream/
SEEKR=~/$run/results/downstream/SEEKR
LPI=~/$run/results/downstream/LPI
RNAfold=~/$run/results/downstream/RNAFold
miranda=~/$run/results/downstream/miranda
LncTar=~/$run/results/downstream/LncTar
filter=/home/cluster/nath/FINAL_SCRIPT/demo/200ntfilter.pl
delnc=/home/cluster/nath/FINAL_SCRIPT/demo/tximport_deseq.R

#Please populate the following variables with their respective local installation (executable) paths

gffcom=/home/cluster/nath/ETENLNC_cluster/tools/gffcompare/gffcompare
cpc=/home/cluster/nath/ETENLNC_cluster/tools/CPC2/bin/CPC2.py
lnctarr=/home/cluster/nath/ETENLNC_cluster/tools/LncTar/LncTar.pl
rnafold=/home/cluster/nath/ETENLNC_cluster/tools/RNAfold
miraanda=/home/cluster/nath/ETENLNC_cluster/tools/miranda
capsule=/home/cluster/nath/ETENLNC_cluster/tools/LPI/
lpic=/home/cluster/nath/ETENLNC_cluster/tools/fasta_data/

echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING SYMLINKS"

ln -s $fastq_query $fastq
ln -s $ref_fa $reffa
ln -s $ref_gtf $refgtf
ln -s $lncRNAs $klncRNA
ln -s $mirna $kmirna
ln -s $proteins $kproteins
ln -s $pc_t $pmrna
touch ~/$run/summary.txt

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
Rscript tximport_deseq.R $base $run varianceStabilizingTransformation

#Fixes left shift
mv -v "$run"\ .csv bak_"$run"\ .csv
echo -n "Transcript_ID", > $run.csv
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
echo -n "Transcript_ID", > $run.csv
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
echo -n "Transcript_ID", > $run.csv
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
echo [`date +"%Y-%m-%d %H:%M:%S"`] "FUNCTIONAL ANNOTATION BY SEEKR COMPLETE"


echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING PROTEIN INTERACTION PREDICTION BY CAPSULE-LPI"
cd $capsule
python3 main.py
cp $lpic/lnc_protein.csv $LPI

cd $LPI
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

echo "ETENLNC RUN SUMMARY" >> ~/$run/summary.txt
echo "PLEASE FIND BELOW THE DETAILS FOR THE PRESENT RUN" >> ~/$run/summary.txt
echo "" >> ~$run/summary.txt
echo "" >> ~/$run/summary.txt
echo "Current Run ID: $run" >> ~/$run/summary.txt
echo "NOVEL LncRNAs" >> ~/$run/summary.txt
echo "Total novel lncRNAs identified in the current run: $p" >> ~/$run/summary.txt
echo "Total novel lncRNAs found to be differentially expressed: $b" >> ~/$run/summary.txt
echo "" >> ~/$run/summary.txt
echo "KNOWN LncRNAs" >> ~/$run/summary.txt
echo "Total known lncRNAs captured in the current run:" $(($n - $p)) "out of" $m >> ~/$run/summary.txt
echo "Total known lncRNAs found to be differentially expressed: $a" >> ~/$run/summary.txt
echo "mRNAs" >> ~/$run/summary.txt
echo "Total mRNAs found to be differentially expressed: $c" >> ~/$run/summary.txt
echo "" >> ~/$run/summary.txt
echo "End_of_run" >> ~/$run/summary.txt

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

