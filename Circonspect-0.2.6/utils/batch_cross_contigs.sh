#! /usr/bin/env bash

# Circonspect, copyright 2010 Florent Angly <florent.angly@gmail.com>, under
# the GNU GPLv3 license. Circonspect generates contig spectra by boostrapping
# assemblies of random subsets of metagenomes.

# batch_cross_contigs.sh is a BASH script to automate creating many cross-contig
# spectra using Circonspect on a multi-core machine or a cluster.

########## CONFIGURATION SECTION ##########
# Edit this section to fit your needs

# Directory where the metagenomes are located
IN_DIR="/home/fangly/Metagenomes"
# FASTA and QUAL file extensions for the metagenomes
FA_SUFFIX=".fa"
QUAL_SUFFIX=".qual"
# Where to write the results
OUT_DIR="/home/fangly/Cross_contigs"
# Path to the Circonspect program
CIRCONSPECT_BIN="/home/fangly/bin/Circonspect"
# Number of random sequences in each sample
SAMPLE_SIZE=10000
# Desired library coverage
COVERAGE=1
# Sequence trimming length
TRIM=100
# Discard sequences below this cutoff
DISCARD=100
# Below, you can specify EITHER A) specific pairs of metagenomes to use, OR B) a
# list of metagenomes to run pairwise. Leave the option you do not want to use
# blank.
# A) Pairs of metagenomes, without file extension (e.g. "metagenome_A metagenome_B")
PAIRS=(
  "metagenome_A metagenome_B"
  "metagenome_A metagenome_C"
)
# B) List of metagenomes, without file extension
LIST=(
)
# for a list of metagenomes, whether or not to generate self vs self cross-contig
# spectra, e.g. metagenome_A vs metagenome_A
INCLUDE_CONTROLS=1

########## END OF CONFIGURATION ##########
# No user-serviceable parts below

# Sanity checks
USAGE="Usage: $0 JOB_ID NUMBER_OF_JOBS\nThis is a batch script to automate generating cross-contig spectra with Circonspect. The computation is split into as many jobs as desired. The different jobs can be run independently, for example on different cores of a CPU or different nodes of a cluster."
NOF_JOBS=$2
if [[ $NOF_JOBS == "" ]]; then echo -e "Error: Missing NUMBER_OF_JOBS argument\n${USAGE}" >&2; exit 0; fi
JOB_NUM=$1
if [[ $JOB_NUM == "" ]]; then echo -e "Error: Missing JOB_ID argument\n${USAGE}" >&2; exit 0; fi
if [[ ! -e $OUT_DIR ]]; then mkdir $OUT_DIR; fi

# Create the list of all pairs to use
NOF_LISTED=${#LIST[*]}
NOF_PAIRS=0
if [[ $NOF_LISTED > 0 ]]; then
  if [[ $INCLUDE_CONTROLS == 1 ]]; then
    OFFSET=0
  else
    OFFSET=1
  fi
  for (( I=1; I<=$NOF_LISTED+$OFFSET; I=$I+1 )); do
    NAME1=${LIST[$I-1]}
    echo "$I: $NAME1"
    for (( J=$I+$OFFSET; J<=$NOF_LISTED; J=$J+1 )); do
      NAME2=${LIST[$J-1]}
      echo "  $J: $NAME2"
      NOF_PAIRS=`expr $NOF_PAIRS + 1`
      PAIRS[$NOF_PAIRS-1]="$NAME1 $NAME2"
    done
  done
fi

# Process all pairs
NOF_PAIRS=${#PAIRS[*]}
for (( I=$JOB_NUM; I<=$NOF_PAIRS; I=$I+$NOF_JOBS )); do
  PAIR=${PAIRS[`expr $I - 1`]}
  NAME1=`echo $PAIR | cut -f 1 -d ' '`
  NAME2=`echo $PAIR | cut -f 2 -d ' '`
  FASTA1=${IN_DIR}/${NAME1}${FA_SUFFIX}
  FASTA2=${IN_DIR}/${NAME2}${FA_SUFFIX}
  QUAL1=${IN_DIR}/${NAME1}${QUAL_SUFFIX}
  QUAL2=${IN_DIR}/${NAME2}${QUAL_SUFFIX}
  BASENAME=${NAME1}_vs_${NAME2}
  OUT=${OUT_DIR}/${BASENAME}.csp
  echo $BASENAME
  if [ -e $QUAL1 ]; then
    QUAL_CMD="-q $QUAL1"
  fi
  if [ -e $QUAL2 ]; then
    QUAL_CMD="$QUAL_CMD $QUAL2"
  fi
  #echo "$CIRCONSPECT_BIN -f $FASTA1 $FASTA2 $QUAL_CMD -k $COVERAGE -s $SAMPLE_SIZE -u $DISCARD -v $TRIM -g -x -e -t -b $BASENAME -d $OUT_DIR > $OUT"
  time $CIRCONSPECT_BIN -f $FASTA1 $FASTA2 $QUAL_CMD -k $COVERAGE -s $SAMPLE_SIZE -u $DISCARD -v $TRIM -g -x -e -t -b $BASENAME -d $OUT_DIR > $OUT
done 

echo "Done!"
exit 1
