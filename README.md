# ETENLNC
ETENLNC (End-To-End-Novel-NonCoding) lncRNA identification and analysis framework

![ETENLNC workflow diagram.](https://github.com/EvolOMICS-TU/ETENLNC/blob/main/ETENLNC.svg?raw=TRUE)

## Installation
ETENLNC can be easily installed and run through a Docker image. Please follow the steps below:
1. Install Docker (https://www.docker.com/products/docker-desktop/)
2. Open up a shell with Docker running in the background (Powershell/cmd for windows, zsh for Mac and terminal for Linux/UNIX)
3. Pull the Docker image using the command:
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

## Running ETENLNC
ETENLNC can be run using the Docker image as follows:
1. Organize all files required for the run under a working directory
2. Open up a shell in the working directory, with Docker running in the background (Powershell/cmd for windows, zsh for Mac and terminal for Linux/UNIX)
3. Run the Docker image using the command:
`docker run -v $(pwd):/user_data/ -w /user_data/ -it prangannath/etenlnc:latest sh`
4. Your working directory has now been mounted under /user_data/. All file paths from now on will be relative to /user_data/ (eg., if you have your ref_genome.fa in the root of the working directory, the path would be modified as /user_data/ref_genome.fa)
5. Run the required ETENLNC docker script using command:
`bash ETENLNC_docker.sh`
6. Enter a suitable RUN ID (run identifier) and follow the instructions to provide your file paths (P.S: The paths must be relative to the docker mounting directory /user_data/)
7. After the run has completed, a directory by your RUN ID name will be generated in your working directory. All results can be found inside this results directory.

## Demo Data
Please find the sample/demo data for ETENLNC here: 
`https://zenodo.org/records/14325721?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVmYTllNWJkLWE4ODUtNGM3OC05YTAxLWM4ZDk5YTljZDZjMCIsImRhdGEiOnt9LCJyYW5kb20iOiI4M2I1YjBlZDQ4MmUyZTIxNDg2YmI0YTFkMWE5MTI1OCJ9.HeB3WrsPduNzMyXjH4x5HfCgmIp4NzYv0P_11XU9lcXC_ZxEaVndP-kD0LDkxwufHVNlroeAhdK33PR51F6jnw`
A detailed guide on running ETENLNC using the demo data can be found in the ETENLNC manual (supplementary to our publication)

## Citation
If you have used ETENLNC for your research, please cite: 
***Nath P, Bhuyan K, Bhattacharyya DK, Barah P. ETENLNC: An end to end lncRNA identification and analysis framework to facilitate construction of known and novel lncRNA regulatory networks. Comput Biol Chem. Published online June 30, 2024. doi:10.1016/j.compbiolchem.2024.108140***
