#!/bin/bash

#SBATCH --run-number=45
#SBATCH --ntasks=50

#####################################################################################################################################################################
#########################################################################     INSTRUCTIONS     ######################################################################
#####################################################################################################################################################################

# 1. Please populate the fastq_query, ref_fa, ref_gtf and lncRNAs with their respective locations
# 2. Please assign the number of threads for the pipeline to use
# 3. Please set a run 'ID' or 'name'. This directory will be created at the current $HOME location
# 4. Please provide valid 'reference_genome' and 'reference_annotation' files for the pipeline. 
# 5. Please provide valid 'lncRNA_fasta' file containing the lncRNAs of the species of interest. 
# 6. Please provide valid 'transcriptome_fast' file for quantification of LncRNAs


#####################################################################################################################################################################
#####################################################################################################################################################################
#######################################################################     START OF SCRIPT     #####################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

zenity --info \
       --title "STARTING ETENLNC" \
       --width 500 \
       --height 100 \
       --text "Welcome to ETENLNC, a one-stop solution for all your lncRNA analysis requirements. Please follow the instructions as provided.
Thank you!"

start_time=`date +%H:%M:%S`
echo "start time:" $start_time

#: << 'END'

#####################################################################################################################################################################
#####################################################################################################################################################################
#############################################################################   PATHS   #############################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

run=`zenity --entry \
       --width 500 \
       --title "Enter Run ID" \
       --text "Please specify an ID for this run. It can be a name/number or experiment/project name."`
       
path="/$run/"

fastq_query=`zenity --file-selection \
	--title "Please select FASTQ files (.fastq.gz)" \
	--filename "/home/${USER}/" \
	--file-filter='*.fastq.gz' \
	--multiple \
	--separator=" "`						#path to your fastq files (fastq.gz)
#echo "\"$fastq_query\" selected."					

ref_fa=`zenity --file-selection \
	--title="Please select Reference Genome file (.fa)" \
	--filename "/home/${USER}/" \
	--file-filter='*.fa'`						#reference_genome (.fa) 
	  
ref_gtf=`zenity --file-selection \
	--title="Please select Reference Annotation file (.gtf)" \
	--filename "/home/${USER}/" \
	--file-filter='*.gtf'` 					#reference_genome (.gtf)   

lncRNAs=`zenity --file-selection \
	--title="Please select known lncRNA fasta file (.fa)" \
	--filename "/home/${USER}/" \
	--file-filter='*.fa'`						#known lncRNA sequences (.fa)  
  
transcripts=`zenity --file-selection \
	--title="Please select Transcriptome file (.fa)" \
	--filename "/home/${USER}/" \
	--file-filter='*.fa'`						#reference transcriptome (.fa)    

deseq=`zenity --file-selection \
	--title="Please select R script file for DEA (.R)" \
	--filename "/home/${USER}/" \
	--file-filter='*.R'`						#R script for Tximport and DESeq2    
 
script=`zenity --file-selection \
	--title="Please select the 'bin' directory created during installation" \
	--filename "/home/${USER}/" \
	--directory`							#bin directory containing downstream analysis tools    

threads=`zenity --entry \
       --width 500 \
       --title "Enter No. of Threads" \
       --text "Please specify the maximum number of processor threads to use for this run."`
       
200nt=`zenity --file-selection \
       --title="Please select 200nt filter perl script (.pl)" \
       --filename "/home/${USER}/" \
       --file-filter='*.pl'`	

#####################################################################################################################################################################
#####################################################################################################################################################################
#################################################################### Start of Pipline workflow_PE ###################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING FASTQ2LNC"

(
# =================================================================
echo "# CREATING DIRECTORIES" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING DIRECTORIES"
mkdir -p ~/$run/results/fastqc/
mkdir -p ~/$run/results/fastp/
mkdir -p ~/$run/results/fastp_reports/
mkdir -p ~/$run/results/featurecounts/
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
mkdir -p ~/$run/results/downstream/
mkdir -p ~/$run/results/downstream/SEEKR
mkdir -p ~/$run/results/downstream/LPI
mkdir -p ~/$run/results/downstream/RNAFold
mkdir -p ~/$run/results/downstream/miranda

mkdir -p ~/$run/query/fastq/
mkdir -p ~/$run/query/ref_fa/
mkdir -p ~/$run/query/ref_gtf/
mkdir -p ~/$run/query/known_lncRNAs/
touch ~/$run/results/lists/list1.txt
touch ~/$run/results/lists/list2.txt

#Permissions
chmod u+x $script

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DIRECTORIES CREATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "ASSIGNING DIRECTORIES"

fastqc=~/$run/results/fastqc/
fastp=~/$run/results/fastp/
fastpr=~/$run/results/fastp_reports/
hisat2=~/$run/results/hisat2/
samtools=~/$run/results/samtools/
fastq=~/$run/query/fastq/
reffa=~/$run/query/ref_fa/
refgtf=~/$run/query/ref_gtf/
index=~/$run/results/index/
fc=~/$run/results/featurecounts/
stringtie=~/$run/results/stringtie/
list1=~/$run/results/lists/list1.txt
list2=~/$run/results/lists/list2.txt
gffcompare=~/$run/results/gffcompare/
gff=~/$run/results/gffcompare/gffcompare-0.11.4.Linux_x86_64/
classes=~/$run/results/lnc_classes/
filters=~/$run/results/filters/
klncRNA=~/$run/query/known_lncRNAs/
blast=~/$run/results/BLAST/
novel=~/$run/results/Predicted_LncRNAs/
base=/home/cluster/nath/quant/
down=~/$run/results/downstream/
SEEKR=~/$run/results/downstream/SEEKR
LPI=~/$run/results/downstream/LPI
RNAfold=~/$run/results/downstream/RNAFold
miranda=~/$run/results/downstream/miranda



#: << 'END'
# =================================================================
echo "2"
echo "# CREATING SYMLINKS" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "CREATING SYMLINKS"

ln -s $fastq_query $fastq
ln -s $ref_fa $reffa
ln -s $ref_gtf $refgtf
ln -s $lncRNAs $klncRNA

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING QUALITY CHECK & PROCESSING"

# =================================================================
echo "5"
echo "# RUNNING FASTP" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING FASTP"

cd $fastq

find $fastq -name "*.fastq.gz" | sort | paste - - | while read A B

do
a=`basename ${A} | sed 's/.sra_1/_1/' | awk -F "." '{print $1}'`
b=`basename ${B} | sed 's/.sra_2/_2/' | awk -F "." '{print $1}'`

# =================================================================
echo "7"
echo "# TRIMMING FILES" ; sleep 2
echo ""
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

# =================================================================
echo "10"
echo "# RUNNING FASTQC" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING FASTQC"

cd $fastp
fastqc *.gz -t $threads
mv *.html "$fastqc"  
mv *.zip "$fastqc" 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QC REPORTS GENERATED"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUALITY CHECK AND PROCESSING DONE"

echo " "
echo " "

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING ALIGNMENT AND ASSEMBLY"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING HISAT2"

# =================================================================
echo "14"
echo "# BUILDING INDICES" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "BUILDING INDICES"

cd $index
hisat2-build -p $threads $reffa/*.fa index


# =================================================================
echo "20"
echo "# RUNNING ALIGNMENT" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING ALIGNMENT"


find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

hisat2 --threads $threads --dta -x $index/index -1 ${A} -2 ${B} -S  $hisat2/$a.sam		
done

#END

# =================================================================
echo "40"
echo "# CONVERTING TO BAM" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "CONVERTING TO BAM"

cd $hisat2

for i in *.sam
do 
echo "Converting $i"
samtools sort -@ "$threads" -o $samtools/$i.bam $i		
echo "$i converted"
done


cd $samtools
#cp *.bam $stringtie

#END


# =================================================================
echo "45"
echo "# RUNNING STRINGTIE" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING STRINGTIE"

#cd $stringtie
for i in `ls *.bam`
do
stringtie -p $threads -G $ref_gtf $i -v -o $i.gtf 
echo "$i.gtf">> $list1
done

stringtie --merge -G $ref_gtf -p $threads -v -o $stringtie/stringtie_merged.gtf $list1  	#--rf-> strandedness, reverse strand; --fr-> forward strand
echo "stringtie_merged.gtf" >> $list2

cp *.bam.gtf $stringtie 
cd $stringtie
cp stringtie_merged.gtf $gffcompare

echo [`date +"%Y-%m-%d %H:%M:%S"`] "ALIGNMENT AND ASSEMBLY DONE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING GENE EXPRESSION ANALYSIS"

cd $samtools
featureCounts -p -T 10 -t gene -g gene_id -a $ref_gtf -o counts.txt -M *.bam

mv *.txt $fc
mv *.summary $fc

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING LNC IDENTIFICATION"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING GFFCOMPARE"

chmod u+x $script/gffcompare/gffcompare
$script/gffcompare/gffcompare -r $refgtf/*.gtf -o $gffcompare/gffannotated.gtf -i $list2
cp gffannotated.gtf.annotated.gtf gffannotated.gtf
cp gffannotated.gtf $classes


# =================================================================
echo "55"
echo "# ISOLATING CLASSES" ; sleep 2
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
cp ~/pipeline/200ntfilter.pl $filters
cp $classes/selected.fa $filters
cp $gffcompare/gffannotated.gtf $filters
cp $classes/ref.fa $filters
perl $200nt 200 selected.fa > gt200.fa

# =================================================================
echo "55"
echo "# RUNNING CODING POTENTIAL ANALYSIS" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING CODING POTENTIAL ANALYSIS"

wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz
gunzip v1.0.1.tar.gz
tar xvf v1.0.1.tar
cd CPC2_standalone-1.0.1 
export CPC_HOME="$PWD"
cd libs/libsvm
gunzip libsvm-3.18.tar.gz 
tar xvf libsvm-3.18.tar
cd libsvm-3.18
make clean && make
cd $filters

./CPC2_standalone-1.0.1/bin/CPC2.py -i gt200.fa -o CPC
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
#bedtools getfasta -fi ref.fa -fo Noncoding200nt.fa -bed Exon_filtered.gtf

#stringtie -eB -G Exon_filtered.gtf -p $threads -v 

#wget https://netix.dl.sourceforge.net/project/rna-cpat/v3.0.0/CPAT-3.0.0.tar.gz
#gunzip CPAT-3.0.0.tar.gz
#tar xvf CPAT-3.0.0.tar

#./CPAT-3.0.0/bin/cpat.py -g Noncoding200nt.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FILTERING DONE"

# =================================================================
echo "58"
echo "# RUNNING BLAST AGAINST KNOWN LNCRNAS" ; sleep 2
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

# =================================================================
echo "80"
echo "# PUTATIVE LNCRNAS IDENTIFIED" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "PUTATIVE LNCRNAS IDENTIFIED"

#: << 'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING QUANTIFICATION OF LNCRNAS"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RUNNING SALMON"

cd $base

echo [`date +"%Y-%m-%d %H:%M:%S"`] "INDEXING THE TRANSCRIPTOME"

salmon index -t $transcripts -i index -k 31

echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFYING...."

find $fastp -name "*_trimmed.fastq.gz" | sort | paste - - | while read A B

do

a=`basename ${A} | awk -F "." '{print $1}' | awk -F "_" '{print $1}'`

salmon quant -i index -l A -1 ${A} -2 ${B} --validateMappings -o $a

done

# =================================================================
echo "85"
echo "# QUANTIFICATION DONE!!!!" ; sleep 2
echo [`date +"%Y-%m-%d %H:%M:%S"`] "QUANTIFICATION DONE!!!!"

: << 'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING TXimport AND DESeq2"

cp $deseq $base

R -f tximport_deseq.R

mv -v res.csv bak_res.csv
echo -n "", > res.csv; cat bak_res.csv >> res.csv #fixes the left shift of column names
rm bak_res.csv

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DE COMPLETE"

END

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FINISHING UP"

cd $novel
p=`cat Novel_LncRNAs.fa | grep '>' | wc -l`
cd $blast
n=`cat $filters/Noncoding200nt.fa | grep '>' | wc -l`
m=`cat $lncRNAs | grep '>' | wc -l`
echo "Putative LncRNAs Identified: $p"
echo "LncRNAs captured from known LncRNAs:" $(($n - $p)) "out of" $m

:<< 'END'

echo [`date +"%Y-%m-%d %H:%M:%S"`] "DOWNSTREAM ANALYSES"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING LNCRNA TARGET PREDICTION BY LncTar"

cd $down

wget http://www.cuilab.cn/lnctarapp/download
unzip download

cd LncTar
path=`pwd`
export PERL5LIB=$path:$PERL5LIB

perl LncTar.pl -p 1 -l $novel/Novel_LncRNAs.fa -m mRNA.txt -d -0.1 -s F -o Targets.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "TARGET PREDICTION COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING FUNCTIONAL ANNOTATION BY SEEKR"

cd $SEEKR
for i in {1..8}
do
seekr_kmer_counts $Novel/Novel_LncRNAs.fa -o kmer_counts_$i.csv -k $i
done

echo [`date +"%Y-%m-%d %H:%M:%S"`] "FUNCTIONAL ANNOTATION BY SEEKR COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING PROTEIN INTERACTION PREDICTION BY CAPSULE-LPI"

python3 main.py
cat example_output.csv | awk -F, '$4==1.0 {print $0}' > new.csv

echo [`date +"%Y-%m-%d %H:%M:%S"`] "PROTEIN INTERACTION PREDICTION BY CAPSULE-LPI COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING RNA FOLD PREDICTION BY RNAFOLD"

cd $script/bin/Capsule_LPI_sources/data/tools
chmod u+x RNAfold
./RNAfold < Novel_LncRNAs.fa
mv *.ps $RNAfold

echo [`date +"%Y-%m-%d %H:%M:%S"`] "RNA FOLD PREDICTION BY RNAFOLD COMPLETE"

echo [`date +"%Y-%m-%d %H:%M:%S"`] "STARTING LNCRNA-miRNA PREDICTION BY miRanda"

cd $script/bin/miranda/
chmod u+x miranda
./miranda mature_high_conf.fa Novel_LncRNAs.fa -out lncRNA-miRNA.txt
mv *.txt $miranda

echo [`date +"%Y-%m-%d %H:%M:%S"`] "LNCRNA-miRNA PREDICTION BY miRanda COMPLETE"

END

echo [`date +"%Y-%m-%d %H:%M:%S"`] "END OF PIPELINE"

# =================================================================
echo "# All finished." ; sleep 2
echo "100"

) |
zenity --progress \
  --title="Progress Status" \
  --text="TEST" \
  --percentage=0 \
  --auto-close \
  --auto-kill

(( $? != 0 )) && zenity --error --text="Error in zenity command."

exit 0

#####################################################################################################################################################################
#########################################################################     END OF SCRIPT    ######################################################################
#####################################################################################################################################################################


















