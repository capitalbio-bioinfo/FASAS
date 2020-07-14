#!/bin/bash

##########################################################################################################################
#
#    Part of FASAS tools script
#
#    Author of pipeEVERY: Ke Zhang
#    Version: 1.1.8
#    For questions, bugs, and suggestions, contact me at ke.zhang@capitalbio.com.
#
##########################################################################################################################

##########################################################################################################################
#    Part of FASAS tools script
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

function GetHelp() {
    echo "
FASAS PDF Organizer
Usage:
    bash $0 -s|--search_path [path]
Options:
    search_path            strings        A FASAS project folder        [none]
    output_path            strings        Output folder                 [none]

"
    exit 0
}

function echoinfo() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][Main]: <Error> $@"
}


test $# -eq 0 && { echoinfo "\033[31mThe parameter is empty, try bash $0 -h\033[0m"; exit 1; }

SEARCHPATH=
OUTPUT=
while true; do
    case "${1:-''}" in
        -s | --search_path) SEARCHPATH=$2; shift 2;;
        -o | --output_path) OUTPUT=$2; shift 2;;
        -h | -help |--help) GetHelp; shift 1;;
        --) shift ; break ;;
        *) break;;
    esac
done

[ ! ${SEARCHPATH} ] && { echoinfo "\033[31mGetOptions: -s | --search_path is necessary!\033[0m"; exit 1; }
[ ! ${OUTPUT} ] && { echoinfo "\033[31mGetOptions: -o | --output_path is necessary!\033[0m"; exit 1; }

test -d ${SEARCHPATH} || { echoinfo "\033[31msearch_path: ${SEARCHPATH}, not found!\033[0m"; exit 1; }
{
    test -d ${OUTPUT} || mkdir -p ${OUTPUT}
} || {
    echoinfo "\033[31mFailed while checking or creating the output directory!\033[0m"
    exit 1
}

#main
for EVERY in `find ${SEARCHPATH} |grep 'pdf$'`; do
    #Get the sample name by the directory name created by FASAS
    name=$(echo $EVERY | perl -ne 'if ($_ =~ /.*\/(.+?)\/step/ or $_ =~ /.*\/(.+?)\/check/){print $1;}')
    #Create new file name
    FILE_BASE_NAME=$(basename $EVERY)
    NEWNAME=$name'_'${FILE_BASE_NAME}
    #Begin copy
    NEWPATH=${OUTPUT}'/'${NEWNAME}
    echo "copy information: cp ${EVERY} -> to -> ${NEWPATH}"
    cp $EVERY ${NEWPATH}
done

echoinfo "Done."
