#!/bin/bash

##########################################################################################################################
#
#    16S Full-Length Amplicon Sequencing Analysis System
#
#    Author of pipeline: Ke Zhang
#    Version: 1.1.8
#    For questions, bugs, and suggestions, contact me at ke.zhang@capitalbio.com.
#
##########################################################################################################################

##########################################################################################################################
#    Base Analysis Script for FASAS
#    Copyright (C) 2019 CapitalBio Corporation
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
##########################################################################################################################

#IDBA_ud has a probability of segmentation error, if IDBA_ud has a segmentation fault, the shell script will rerun it
#In order to avoid perl receiving an error signal directly when calling IDBA_ud, package IDBA_ud with a shell script.

FastaFile=''
OutputFolder=''

function GetHelp {
    echo "
Full-length 16S Assemble Program
Usage:
    bash $0 -f|--fasta_file [input file] -o|--output_folder [output folder] -p|--program [idba_ud || cap3]
Parameter
    FastaFile       path    input file, FASTA format. [required]
    OutputFolder    path    output folder. [required]
    AssembleProgram strings program used for assembly. [required]
"
    exit 0;
}

function echowarn() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][Main]: <Warn> $@"
}

function echoinfo() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][Main]: <Info> $@"
}

function echoerr() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][Main]: <Error> $@"
}

##if opt == 0
test $# -lt 6 && { echo -e "\033[31mThe program lacks the necessary parameters\033[0m"; GetHelp; exit 1; }

while true; do
    case "${1:-''}" in
        -f | --fasta_file) FastaFile=$2; shift 2;;
        -o | --output_folder) OutputFolder=$2; shift 2;;
        -p | --program) AssembleProgam=$2; shift 2;;
        -h | -help |--help) GetHelp; shift 1;;
        --) shift ; break ;;
        *) break;;
    esac
done

test ${#FastaFile} -eq 0 && { echo "FastaFile must have a value"; exit 1; }
test -r ${FastaFile} || { echo "ERROR: ${FastaFile} cannot be read!"; exit 1; }

function RunCap3 {
    CapFile=${FastaFile%.*}
    CapReuslt=${CapFile}'_r1.cap.out'
    CapFile=${CapFile}'_r1.fa'
    {
        cat ${FastaFile} | perl -ne '$seq=<>;$head=<>;$seq2=<>;print "$_$seq"' > ${CapFile}
    } || {
        echo "ERROR: Failed to Assemble $FastaFile, When Creating a Cap3 Sequence File"
        exit 1
    }
    cap3 ${CapFile} > ${CapReuslt}
}

function RunIDBA {
    #run idba_ud, and get log file
    #But, idba_ud is calculating the possibility of aborting the program due to a segmentation error
    #Note: the big probability of this error occurs before the scaffold is generated
        #small probability before the log file is generated
    #So. I added a 3 loop to the code below to rerun idba_ud
    idba_ud -r ${FastaFile} -o ${OutputFolder} --mink 26 --maxk 46 --step 4 --num_threads 1 &> /dev/null
    #check log file
    ReRunlimit=1
    LogFile=${OutputFolder}'/log'
    if [ ! -e ${LogFile} ]; then
        echowarn "${FastaFile}, not found LogFile"
        while [ 1 ]; do
            idba_ud -r ${FastaFile} -o ${OutputFolder} --mink 26 --maxk 46 --step 4 --num_threads 1 &> /dev/null
            if [[ -e ${LogFile} ]]; then
                echoinfo "${FastaFile}, jump out of the loop in normal state"
                break
            elif [[ ${ReRunLimit} -gt 3 ]]; then
                echowarn "${FastaFile}, jump out of the loop because of the number os loops"
                break
            fi
            let ReRunLimit=${ReRunLimit}+1
        done
    fi
}

#main
case "${AssembleProgam}" in
    cap3) RunCap3;;
    idba_ud) RunIDBA;;
esac

