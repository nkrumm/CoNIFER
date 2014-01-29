#!/bin/bash
# Trap interrupts and exit instead of continuing the loop
trap "echo Exited!; exit;" SIGINT SIGTERM

COMPONENTS_REMOVED=70
CONIFER_ANALYSIS_FILE="ESP.UW.SVD${COMPONENTS_REMOVED}.QC.hdf5"
QC_REPORT="ESP.UW.QC.SVD${COMPONENTS_REMOVED}.txt"

PROBES="/net/eichler/vol8/home/nkrumm/EXOMES/ESP2000/UW.probes.txt"
#SAMPLES="/net/eichler/vol8/home/nkrumm/EXOMES/ESP2000/UW.samples.txt"
SAMPLES="/net/eichler/vol8/home/nkrumm/EXOMES/ESP2000/Calling/UW/SVD70/UW.samples.passQC.txt"
LOGLEVEL="INFO"


CHRS=`seq 1 24`
FILES=""

TMPDIR=/var/tmp/`whoami`
mkdir -p $TMPDIR
EXCEPTION=0
for CHR in $CHRS; 
do
	echo $CHR
	CHR_FILE="${TMPDIR}/chr${CHR}.hdf5"
	QC_report="${TMPDIR}/chr${CHR}.qc.txt"

	python 02_create_conifer_file.py -O $CHR_FILE --components_removed $COMPONENTS_REMOVED --chromosomes $CHR --probes $PROBES --samples $SAMPLES --loglevel $LOGLEVEL --QC_report=$QC_report
	# chceck return value, exceptions will return 1, and we exit the loop
	if [ $? == 1 ]; then
        echo "Exception Caught!"
        EXCEPTION=1
        break
    else
    	#otherwise file was successfully made, append it to list of files for merging later
        FILES="${FILES} ${CHR_FILE}"
    fi
done

if [ $EXCEPTION != 1 ]; then
	echo "Merging all files..."
	python 03_merge_conifer_files.py --outfile $CONIFER_ANALYSIS_FILE --infiles $FILES
	echo "Making concatenated QC file"
	cat ${TMPDIR}/chr*.qc.txt > $QC_REPORT
fi
echo "Cleaning up..."
rm -f $FILES
rm -f ${TMPDIR}/chr*.qc.txt
echo "Done."