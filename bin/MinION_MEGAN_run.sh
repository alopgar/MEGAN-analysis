#!/bin/bash
#SBATCH -N 1      #solicita un modulo
#SBATCH -n 1      #numero de tareas
#SBATCH -c 24     #Numero de procesadores por tarea
#SBATCH --partition=thinnodes
#SBATCH -t 09:00:00      #tiempo de ejecucion hh:mm:ss
#SBATCH -J DIAMOND     #nombre del proceso

## No indicamos sbatch --mem=40GB porque estamos pidiendo un nodo de thinnodes,cola-corta en exclusiva (120GB). Dejamos que use toda la memoria.
## Default time: 10h

echo -e "DIAMOND + MEGAN pipeline for Nanopore MinION metagenome sequencing.
Project: $WDIR
Analyzing sample "$(awk '/'$1'/{ print NR; }' $IDfile)" of "$(< $IDfile wc -l)".
Input file: $1.fastq"

####################################################################
# PRINSEQ for quality & length trimming
####################################################################
module load perl
perl $PRINSEQDIR/prinseq-lite.pl -fastq $WDIR/1_Rawdata/$1.fastq -min_qual_mean 7 -min_len 300 \
	-out_good $WDIR/1_Rawdata/$1"_trim" -out_bad $WDIR/1_Rawdata/trim_bad/$1"_bad" -log $WDIR/1_Rawdata/trim_logs/$1.log

####################################################################
# DIAMOND for .fastq to .daa converting 
####################################################################
module load gcc/6.4.0 diamond/0.9.22

##OPTIONS:
## blastx is the alignment option
## -query points to the file to convert
## -db is database (instructions to get this ask Bea)
## -daa output file
## -F 15 (to consider frameshift errors from Nanopore reads)
## -range-culling --top 10 (consider only alignments with a bitscore in the best 10% bit score)

diamond blastx -F 15 --range-culling --top 10 --query $WDIR/1_Rawdata/$1"_trim.fastq" --db $DMND --daa $WDIR/2_files.daa/$1"_trim.daa"

####################################################################
# MEGAN for .daa to .rma converting
####################################################################
##OPTIONS
## -i				input_file
## -o				ouput file
## -a2t				Map alignment against accesion map of the *proteins* for *taxonomy classification* (Mapping files for current NCBI-nr protein database (not containing GI numbers). Protein to taxonomy http://ab.inf.uni-tuebingen.de/data/software/megan6/download/welcome.html)
## -a2interpro2go	Alignment from *proteins* to *functional information* (Mapping files for current NCBI-nr protein database (not containing GI numbers). Protein to InterPro http://ab.inf.uni-tuebingen.de/data/software/megan6/download/welcome.html)
## -lg				Use long reads
## -alg				Algorithm to use for mapping. Use longReads specific for nanopore
## -ram				Set the read assignment mode. Default value: readCount. Legal values: readCount, readLength, alignedBases, readMagnitude
## -ms				Min score. Default value: 50.0.

$MEGANDIR/tools/daa2rma -i $WDIR/2_files.daa/$1"_trim.daa" -o $WDIR/3_files.rma/$1"_trim.rma" -lg -alg longReads -a2t $MEGAN_a2t -a2interpro2go $MEGAN_ip2g

####################################################################
# MEGAN for taking taxonomic & functional info from .rma
####################################################################
#Converts rma file into a readable file with taxonomy and number of reads per taxon.
$MEGANDIR/tools/rma2info -i $WDIR/3_files.rma/$1"_trim.rma" -c2c Taxonomy -p -r -mro > $WDIR/4_files.rmainfo/$1"_trim.rma.RA.txt"

#Converts rma file into a readable file with Interpro2go functional info per read (2nd line is for taxonomy of each read)
$MEGANDIR/tools/rma2info -i $WDIR/3_files.rma/$1"_trim.rma" -r2c INTERPRO2GO -p > $WDIR/4_files.rmainfo/$1"_trim.rma.IP2G.txt"
$MEGANDIR/tools/rma2info -i $WDIR/3_files.rma/$1"_trim.rma" -r2c Taxonomy -p -r -mro > $WDIR/4_files.rmainfo/$1"_trim.rma.readtotax.txt"
