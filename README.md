# ETENLNC
ETENLNC (End-To-End-Novel-NonCoding) lncRNA identification and analysis framework

## Installation
ETENLNC can be easily installed and run through a Docker image. Please follow the steps below:
1. Install Docker (https://www.docker.com/products/docker-desktop/)
2. Open up a shell with Docker running in the background (Powershell/cmd for windows, zsh for Mac and terminal for Linux/UNIX)
3. Type the following:
   `docker pull prangannath/etenlnc:latest`
4. Sit back, relax. Let it download.
ETENLNC has been successfully installed in your system. Please refer to the next section for instructions on running ETENLNC.

## Input Files
The following inputs are required for running ETENLNC:
* F1. Paired-ended FASTQ files
* F2. A reference-genome for the concerned species (.fasta or .fa format)
* F3. A reference annotation containing all genes and transcripts (.gtf or .gff3 format)
* F4. Known lncRNA sequences (.fasta or .fa format)
* F5. Known protein-coding transcripts (mRNA) sequences (.fasta or .fa format)
* F6. A metadata file having sample-condition information (Samples.txt)

Additional Files for downstream analyses:
* F7. miRNA sequences (.fasta or .fa format)
* F8. protein sequences (.fasta or .fa format)

## Scripts (for Docker)
* ETENLNC_docker.sh: Runs entire ETENLNC pipeline and all additional downstream modules. Requires all files (F1-F8)
* ETENLNC_docker_DEO.sh: Runs only quantification and differential expression of known lncRNAs and mRNAs. Requires files F1, F4, F5 and F6
* ETENLNC_docker_DO.sh: Runs quantification, differential expression of known lncRNAs, mRNAs and downstream predictions to construct networks. Requires files F1, F4-F8
* ETENLNC_docker_IDE.sh: Runs identification of novel lncRNAs, quantification and differential expression of known and novel lncRNAs, mRNAs. Requires files F1-F6
* ETENLNC_docker_IO.sh: Runs identification of novel lncRNAs alone. Requires files F1-F4

## Citation
If you have used ETENLNC for your research, please cite: 
***Nath P, Bhuyan K, Bhattacharyya DK, Barah P. ETENLNC: An end to end lncRNA identification and analysis framework to facilitate construction of known and novel lncRNA regulatory networks. Comput Biol Chem. Published online June 30, 2024. doi:10.1016/j.compbiolchem.2024.108140***
