#!/bin/bash

#Setting the working directory
wd=`pwd`

#Installing all aptitude, pip and R dependencies
apt-get update && apt-get upgrade -y
for i in python3 python3-pip apt-utils wget make fastp hisat2 fastqc stringtie zenity samtools gffread subread salmon unzip r-base ncbi-blast+ r-cran-biocmanager r-cran-stringr r-bioc-deseq2 r-bioc-tximportData r-cran-ggplot2 r-cran-dplyr r-cran-readr r-bioc-tximport r-cran-gplots seqtk r-cran-rcurl libmariadb-dev gcc python3-dev libcogl-pango-dev libcairo2-dev libtool linux-headers-amd64 musl-dev libffi-dev libssl-dev libjpeg-dev zlib1g-dev seekr; do apt-get install -y $i; done
python3 -m pip install --upgrade pip
python3 -m pip install six torch numpy biopython pandas Cython tqdm seekr --no-cache-dir
Rscript -e 'BiocManager::install(c("tximportData", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment", "genefilter", "geneplotter", "DESeq2")'
Rscript -e 'install.packages("tidyverse")'

#Installing tools (Please make sure you have the tools directory from our Github)
chmod -R 777 tools/LPI/data/tools
cd /tools/CPC2/libs/libsvm/libsvm-3.18/
make clean && make
cp tools/Genecare1/Genecare /etc/perl/

#Printing out the installation paths to be provided during running 

echo "Please save the following installation paths. These have to be provided during running the ETENLNC framework"

echo "
gffcom = $pwd/tools/gffcompare/gffcompare
cpc = $pwd/tools/CPC2/bin/CPC2.py
lnctarr = $pwd/tools/LncTar/LncTar.pl
rnafold = $pwd/tools/RNAfold
miraanda = $pwd/tools/miranda
capsule = $pwd/tools/LPI/
lpic = $pwd/tools/fasta_data/"

echo "Please save the following script paths. These have to be provided during running the ETENLNC framework"

echo "
filter = $pwd/scripts/200ntfilter.pl
delnc = $pwd/scripts/tximport_deseq.R"

echo "ETENLNC installation complete!"

exit

