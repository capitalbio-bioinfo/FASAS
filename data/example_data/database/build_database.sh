#!/bin/bash
if [ $(which bowtie2-build 2>/dev/null) ] ; then
    echo "Found bowtie2-build command."
else
    echo "bowtie2-build command not found!"
    echo "exit."
    exit
fi

if [ $(which makeblastdb 2>/dev/null) ] ; then
    echo "Found makeblastdb command."
else
    echo "makeblastdb command not found!"
    echo "exit"
    exit
fi

if [ $(which makembindex 2>/dev/null) ] ; then
    echo "Found makembindex command."
else
    echo "makembindex command not found!"
    echo "exit"
    exit
fi

if [ $1 -eq '' ] ; then
    echo "Requires a parameter, INPUT_FILE."
    echo "exit."
    exit
fi

test -e $1 || { echo "INPUT_FILE: ${INPUT_FILE} is not found!"; echo "exit."; exit; }

INPUT_FILE=$1
INPUT_FOLDER=$(dirname ${INPUT_FILE})
INPUT_BASENAME=$(basename ${INPUT_FILE})
cd ${INPUT_FOLDER}

echo "Build bowtie2 index."
bowtie2-build ${INPUT_BASENAME} ${INPUT_BASENAME} >/dev/null
echo "Build blastn index."
makeblastdb -in ${INPUT_BASENAME} -input_type fasta -dbtype nucl
makembindex -iformat blastdb -input ${INPUT_BASENAME}

echo "Done."
