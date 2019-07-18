#!/bin/bash
#SBATCH -t 01:30:00

usage(){
cat << _EUs_
$(basename "$0") [PATHS]... [OPTIONS]... -- Initiates a pipeline for Nanopore MinION metagenome analysis using DIAMOND + MEGAN.\n
\nDESCRIPTION:\n
	This script performs DIAMOND + MEGAN analysis of .fastq files from .fast5 Nanopore MinION data. It launches multiple tasks to remote server queue,\n
	one per input sample. Modifications on default parameters (trimming software, DIAMOND / MEGAN params...) must be done INSIDE this code (see script\n
	MinION_MEGAN_run.sh).\n
\nPATHS:\n
	First mandatory argument -> Input Path. Where .fastq or .fastq.gz files are stored.\n
	Second mandatory argument -> Output Path. Working directory. Output folders and results will be stored here.\n
\nPARAMETERS:\n
	-i, --init:\n
	\tThis flag indicates that initial files and configs must be prepared in output dir. If NOT SPECIFIED, NO FASTQ SYMLINKS nor ID FILE CREATION\n
	\twill be done.\n
	-F, --final:\n
	\tThis flag enables the use of .log files for creating final stats files:\n
	\t\t-Trimming stats (Prinseq-lite)\n
	\t\t-Running time for DIAMOND files\n
	\t\t-Running time for MEGAN files\n
	\tIf this flag is indicated, DIAMOND + MEGAN processing will NOT be done.\n
	-t, --tempdir:\n
	\tCustomize temp directory where input data will be copied and unzipped.If NOT SPECIFIED, input data will be copied in [Input_path]/temp.\n 
	-h, --help:\n
	\tPrints this screen.\n
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

export PTH=$1
export PTH_data=$PTH/temp
export WDIR=$2; shift 2
ini=
final=0

OPTS=`getopt -o iFt:h --long init,final,tempdir:,help -- "$@"`
eval set -- "$OPTS"

while true; do
	case $1 in
		-i | --init)
			ini=1; shift ;;
		-F | --final)
			final=1; shift ;;
		-t | --tempdir)
			export PTH_data=$2; shift 2 ;;
		-h | --help)
			echo -e $(usage) | less ; exit ;;
		--) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
	esac
done

if [ ! -d $WDIR/1_Rawdata ];then mkdir -p $WDIR/1_Rawdata; fi

echo "------ INPUT DIRECTORY: "$PTH
echo "------ TEMP DIRECTORY: "$PTH_data
echo "------ OUTPUT DIRECTORY: "$WDIR

####################################################################
# RUNNING SCRIPT
####################################################################

## START: IDs FILE & LINKS TO INPUT FILES (if parameter 'ini' specified)
export IDfile=$WDIR/INP_files.txt

if [ "$ini" == 1 ]; then
	echo "INI"
	if [ ! -d $PTH_data ];then mkdir -p $PTH_data; fi
	if [ -f $WDIR/INP_files_all.txt ]; then mv $WDIR/INP_files_all.txt $WDIR/INP_files_done.txt
	else touch $WDIR/INP_files_done.txt; fi
	ls $PTH/*fastq* | xargs -n1 basename | awk '{gsub(".fastq",""); gsub(".gz",""); print }' > $WDIR/INP_files_all.txt
	diff $WDIR/INP_files_all.txt $WDIR/INP_files_done.txt | grep '<' | sed 's/<\ //' > $WDIR/INP_files.txt
	export Input=$(awk '{ print $1 }' $IDfile)
	for f in $Input; do
		echo "Copying file "$f" from INPUT DIR to TEMP DIR"
		cp $PTH/$f* $PTH_data
		if [ -f $PTH_data/$f.fastq.gz ]; 
			echo "File .gz --- gunzipping..."
			then gunzip $PTH_data/$f.fastq.gz
		fi
		ln -s $PTH_data/$f.fastq $WDIR/1_Rawdata
		echo "Copying done --- Symlink created in 1_Rawdata"
	done
else export Input=$(awk '{ print $1 }' $IDfile); fi

## TRIMMING FOLDERS:
if [ ! -d $WDIR/1_Rawdata/trim_bad ];then mkdir $WDIR/1_Rawdata/trim_bad; fi
if [ ! -d $WDIR/1_Rawdata/trim_logs ];then mkdir $WDIR/1_Rawdata/trim_logs; fi

## DIAMOND + MEGAN FOLDERS:
if [ ! -d $WDIR/2_files.daa ];then mkdir $WDIR/2_files.daa; fi
if [ ! -d $WDIR/3_files.rma ];then mkdir $WDIR/3_files.rma; fi
if [ ! -d $WDIR/4_files.rmainfo ];then mkdir $WDIR/4_files.rmainfo; fi

## EXE:
if [ "$final" == 0 ]; then
	for file in $Input; do
		sbatch ./bin/MinION_MEGAN/MinION_MEGAN_run.sh $file
	done
fi

## FINALSTATS:
if [ "$final" == 1 ]; then
	sh ./bin/MinION_MEGAN/MinION_MEGAN_finalstats.sh
fi
