#!/bin/bash
#SBATCH -t 01:30:00

usage(){
cat << _EUs_
$(basename "$0") [OPTIONS]... -- Initiates a pipeline for Nanopore MinION metagenome analysis using DIAMOND + MEGAN.\n
\nDESCRIPTION:\n
	\tThis script performs DIAMOND + MEGAN analysis of .fastq files from .fast5 Nanopore MinION data. It launches multiple tasks to remote \n
	\tserver queue, one per input sample. Modification of input/output paths must be done inside this code, although implementation of new \n
	\tcommand arguments is being developed. Modifications on default parameters (trimming software, DIAMOND/MEGAN params...) must be done \n
	\tINSIDE this code (see script MinION_MEGAN_run.sh).\n
\nPARAMETERS:\n
	\t-i, --init:\n
		\t\tThis flag indicates that initial files and configs must be prepared in output dir.\n 
		\t\tIf NOT SPECIFIED, NO FASTQ SYMLINKS nor ID FILE CREATION will be done.\n
	\t-F, --final:\n
		\t\tThis flag enables the use of .log files for creating final stats files:\n
			\t\t\t-Trimming stats (Prinseq-lite)\n
			\t\t\t-Running time for DIAMOND files\n
			\t\t\t-Running time for MEGAN files\n
		\t\tIf this flag is indicated, DIAMOND + MEGAN processing will NOT be done.\n
	\t-h, --help:\n
		\t\tPrints this screen.\n
_EUs_
}

####################################################################
# INITIAL VARIABLES
####################################################################
export PRINSEQDIR=$STORE/ENVS/prinseq-lite-0.20.4
export MEGANDIR=$STORE/ENVS/megan

export DMND=$LUSTRE/NCBI_db/nr.dmnd
export MEGAN_a2t=$LUSTRE/NCBI_db/prot_acc2tax-May2017.abin
export MEGAN_ip2g=$LUSTRE/NCBI_db/acc2interpro-Nov2016XX.abin

export PTH=/mnt/netapp2/Store_uni/home/otras/ini/ogr/SHARED_DATA/METALGENmetagenomas
export PTH_data=$STORE/Data/FULL_METAGENOME/METALGEN_MINION
export WDIR=$LUSTRE/Results/Metagenome/METALGEN
ini=
final=0

OPTS=`getopt -o i::F::h --long init::,final::,help -- "$@"`
eval set -- "$OPTS"

while true; do
	case $1 in
		-i | --init)
			ini=1; shift 2 ;;
		-F | --final)
			final=1; shift 2 ;;
		-h | --help)
			echo -e $(usage) | less ; exit ;;
		--) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
	esac
done

if [ ! -d $PTH_data ];then mkdir -p $PTH_data; fi
if [ ! -d $WDIR/1_Rawdata ];then mkdir -p $WDIR/1_Rawdata; fi

####################################################################
# RUNNING SCRIPT
####################################################################

## START: IDs FILE & LINKS TO INPUT FILES (if parameter 'ini' specified)
export IDfile=$WDIR/INP_files.txt

if [ "$ini" == 1 ]; then
	echo "INI"
	if [ -f $WDIR/INP_files_all.txt ]; then mv $WDIR/INP_files_all.txt $WDIR/INP_files_done.txt
	else touch $WDIR/INP_files_done.txt; fi
	ls $PTH/*fastq* | xargs -n1 basename | awk '{gsub(".fastq",""); gsub(".gz",""); print }' > $WDIR/INP_files_all.txt
	diff $WDIR/INP_files_all.txt $WDIR/INP_files_done.txt | grep '<' | sed 's/<\ //' > $WDIR/INP_files.txt
	export Input=$(awk '{ print $1 }' $IDfile)
	for f in $Input; do
		cp $PTH/$f* $PTH_data
		if [ -f $PTH_data/$f.fastq.gz ]; then gunzip $PTH_data/$f.fastq.gz; fi
		ln -s $PTH_data/$f.fastq $WDIR/1_Rawdata
	done
else export Input=$(awk '{ print $1 }' $IDfile); fi

## TRIMMING FOLDERS:
if [ ! -d $WDIR/1_Rawdata/trim_bad ];then mkdir $WDIR/1_Rawdata/trim_bad; fi
if [ ! -d $WDIR/1_Rawdata/trim_logs ];then mkdir $WDIR/1_Rawdata/trim_logs; fi

## EXE:
if [ "$final" == 0 ]; then
	for file in $Input; do
		sbatch MinION_MEGAN_run.sh $file
	done
fi

## FINALSTATS:
if [ "$final" == 1 ]; then
	sbatch MinION_MEGAN_finalstats.sh
fi
