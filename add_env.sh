#!/bin/bash

##########################################################################################################################
#   File name: add_env.sh
#   Auther: Chen
#   Version: 1.1.8
#   Data: 2019.11.29
#   Function: Add FASASHome and PERL5LIB to the shell environment
##########################################################################################################################

##########################################################################################################################
#    add_env.sh
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

#This script needs to use bash
#Here, use the contents of the proc folder in the Linux system to determine whether the bash is currently used.
QScriptEXE=$(basename `readlink /proc/$$/exe`)
( test ${QScriptEXE} = 'bash' && echo "[$(date +'%Y%m%d-%H:%M:%S')][add_env]: <Info> Sure, this is bash." ) || { echo -e "[$(date +'%Y%m%d-%H:%M:%S')][add_env]: <Error> \033[31mError: It is recommended to use bash\033[0m"; exit 1; }

#define some functions that output information
function echowarn() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][add_env]: <Warn> $@"
}

function echoinfo() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][add_env]: <Info> $@"
}

function echoerr() {
    echo -e "[$(date +'%Y%m%d-%H:%M:%S')][add_env]: <Error> $@"
}

sleep 1
MYPath=$(dirname `readlink -f ${BASH_SOURCE[0]}`)
MYProgram=${MYPath}'/run_analysis.sh'
if [ ! -e ${MYProgram} ]
then
    echoerr "No 16S-FASAS script found in ${MYPath}."
    echoinfo "Exit."
    exit 1
fi

sleep 1
BashrcPATH=${HOME}'/.bashrc'
FASASPERL5LIB=${MYPath}'/perl'

#export
export FASASHome=${MYPath}
export PERL5LIB=${FASASPERL5LIB}:$PERL5LIB

sleep 1
#check bashrc and add env variable to ~/.bashrc
if [ -s ${BashrcPATH} ]
then
    echoinfo "Find ~/.bashrc and size of ~/.bashrc gt zero."
    if [[ $(cat ${BashrcPATH}|perl -ne 'BEGIN{undef $/} /FASASHome/ ? print 1 : print 0') == 1 ]]
    then
        echoinfo "There are FASASHome environment variable in ${BashrcPATH}"
    else
        echoinfo "Prepare to append environment variable to ~/.bashrc"
        echo -en "\n##16S-FASAS environment variable\n" >> ${BashrcPATH}
        echo -en "export FASASHome=${MYPath}\n" >> ${BashrcPATH}
        echo -en "export PERL5LIB=${PERL5LIBPATH}:\$PERL5LIB\n" >> ${BashrcPATH}
    fi
else
    echoinfo "~/.bashrc not found or size of ~/.bashrc eq 0."
    echoinfo "Create ~/.bashrc and write environment variable to ~/.bashrc"
    echo "Do you approve create ~/.bashrc ?[yes|no]"
    read CreateSignal
    while [[ (${CreateSignal} != "yes") && (${CreateSignal} != "Yes") && (${CreateSignal} != "YES") &&
        (${CreateSignal} != "no") && (${CreateSignal} != "No") && (${CreateSignal} != "NO") ]]
    do
        echo "Please answer 'yes' or 'no':"
        read CreateSignal
    done
    if [[ $(echo -n ${CreateSignal}|perl -ne '/yes/i ? print 1 : print 0') == 1 ]]
    then
        echo -en "\n##16S-FASAS environment variable\n" >> ${BashrcPATH}
        echo -en "export FASASHome=${MYPath}\n" >> ${BashrcPATH}
        echo -en "export PERL5LIB=${PERL5LIBPATH}:\$PERL5LIB\n" >> ${BashrcPATH}
    else
        echowarn "~/.bashrc no create."
        echowarn "Exit."
        exit 2
    fi
fi
