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
#    create_conda_env.sh
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

CondaLink='https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh'

CondaPATH=$(which conda 2> /dev/null)
CondaStatus=$?

if [[ ${CondaStatus} -eq 0 ]]
then
    echo "Conda has been installed, skip"
else
    echo "Conda not found! Prepear download and install Miniconda2-latest-Linux-x86_64!"
    echo "Do you want to download and install miniconda? [yes|no]"
    read CondaSignal
    while [[ (${CondaSignal} != "yes") && (${CondaSignal} != "Yes") && (${CondaSignal} != "YES") &&
        (${CondaSignal} != "no") && (${CondaSignal} != "No") && (${CondaSignal} != "NO") ]]
    do
        echo "Please answer 'yes' or 'no':"
        read CondaSignal
    done
    if [[ $(echo -n ${CondaSignal}|perl -ne '/yes/i ? print 1 : print 0') != 1 ]]
    then
        echo "The operation was rejected by the user, exit."
        exit 2
    fi

    echo "Download Miniconda2-latest-Linux-x86_64"
    DownSitePATH=${HOME}'/Miniconda2-latest-Linux-x86_64.sh'
    DownTools=('axel' 'curl' 'wget')
    for DownToolIndex in {0..2}
    do
        UseTool=$(which ${DownTools[${DownToolIndex}]} 2> /dev/null)
        if [[ ! ${UseTool} =~ ^$ ]]
        then
            break
        fi
    done
    if [[ ${UseTool} =~ ^$ ]]
    then
        echo "No download tools available in the environment: wget, curl or axel."
        echo "exit."
        exit
    elif [[ ${UseTool} =~ wget ]]
    then
        {
            wget -O ${DownSitePATH} ${CondaLink}
        } || {
            echo "ERROR: Return $? when downloading the Miniconda installation script using wget."
            echo "exit."
            exit 1
        }
    elif [[ $UseTool =~ curl ]]
    then
        {
            curl ${CondaLink} > ${DownSitePATH}
        } || {
            echo "ERROR: Return $? when downloading the Miniconda installation script using curl."
            echo "exit."
            exit 1
        }
    elif [[ ${UseTool} =~ axel ]]
    then
        {
            axel -o ${DownSitePATH} ${CondaLink}
        } || {
            echo "ERROR: Return $? when downloading the Miniconda installation script using axel."
            echo "exit."
            exit 1
        }
    fi
    echo "Start installing Miniconda."
    bash ${DownSitePATH}
    if [ $? -ne 0 ]
    then
        echo "ERROR: Return $? when installing Miniconda"
        echo "exit"
        exit 1
    fi
    #check conda install
    echo "Please enter the path to the conda program"
    echo "Example: /home/user/miniconda2/bin/conda"
    read CondaInstallPATH
    while [[ $(echo ${CondaInstallPATH}|perl -ne '/conda/ ? print 1 : print 0;') != 1 ]]
    do
        echo "Please enter the path to the conda program"
        echo "Just like: /home/user/miniconda2/bin/conda"
        read CondaInstallPATH
    done
    if [ ! -x ${CondaInstallPATH} ]
    then
        echo "${CondaInstallPATH} does not exist or is not an executable program"
        echo "exit"
        exit 1
    fi
    CondaPATH=${CondaInstallPATH}
fi

CondaVersion=$(${CondaPATH} --version 2>&1 | perl -ne 'BEGIN{undef $/} /(\d+\.\d+)/; $1 < 4.7 ? print 1 : print 0')
if [[ ${CondaVersion} -eq 1 ]]
then
    echo "Conda has beem installed, but version lt 4.7"
    echo "Do you want to update conda? optional [yes|no]"
    read CondaUpdate
    while [[ (${CondaSignal} != "yes") && (${CondaSignal} != "Yes") && (${CondaSignal} != "YES") &&
        (${CondaSignal} != "no") && (${CondaSignal} != "No") && (${CondaSignal} != "NO") ]]
    do
        echo "Please answer 'yes' or 'no':"
        read CondaUpdate
    done
    if [[ $(echo -n ${CondaUpdate}|perl -ne '/yes/i ? print 1 : print 0') != 1 ]]
    then
        echo "Version of Conda is important!"
        echo "In order to avoid accidents, it is recommended that the version of conda is above 4.7!!!"
        echo "exit."
        exit 2
    fi
    {
        conda update conda
    } || {
        echo "ERROR: Return $? when updating conda"
        echo "exit."
        exit 1
    }
fi

#conda install env
MyName=$(whoami)
CondaChown=$(stat -c %U ${CondaPATH})
if [ ${MyName} == ${CondaChown} ]
then
    echo "Found the conda command"
    echo "Do you want to use conda's default software sources, such as conda-forge and bioconda? [yes|no]"
    echo "Note: Conda's default software source server is in the United States."
    read CondaConfig
    while [[ (${CondaConfig} != "yes") && (${CondaConfig} != "Yes") && (${CondaConfig} != "YES") && 
        (${CondaConfig} != "no") && (${CondaConfig} != "No") && (${CondaConfig} != "NO") ]]
    do
        echo "Please answer 'yes' or 'no':"
        read CondaConfig
    done
    if [[ $(echo -n ${CondaConfig} | perl -ne '/yes/i ? print 1 : print 0') == 1 ]]
    then
        echo "Add r, bioconda and conda-forge to the conda configuration"
        ${CondaPATH} config --add channels r
        ${CondaPATH} config --add channels bioconda
        ${CondaPATH} config --add channels conda-forge
    else
        echo "Using current Conda configuration."
    fi

    {
        echo "Execute: conda create -n FASAS command ......"
        ${CondaPATH} create --yes -n FASAS trimmomatic cutadapt=2.3 r-ggplot2 pigz=2.3.4 bowtie2=2.3.5 r-optparse r-pheatmap perl=5.26.2 perl-file-which perl-parallel-forkmanager
        echo "Execute: conda install --no-deps -n FASAS blast ......"
        ${CondaPATH} install --yes --no-deps -n FASAS blast
    } || {
        echo "An error occurred during the execution of the conda command."
        echo "If you encounter a network connection type error, please try a few more times or modify the configuration of conda."
        exit
    }
else
    echo "Found the conda command"
    echo "Do you want to use conda's default software sources, such as conda-forge and bioconda? [yes|no]"
    echo "Note: Conda's default software source server is in the United States."
    read CondaConfig
    while [[ (${CondaConfig} != "yes") && (${CondaConfig} != "Yes") && (${CondaConfig} != "YES") &&
        (${CondaConfig} != "no") && (${CondaConfig} != "No") && (${CondaConfig} != "NO") ]]
    do
        echo "Please answer 'yes' or 'no':"
        read CondaConfig
    done
    if [[ $(echo -n ${CondaConfig} | perl -ne '/yes/i ? print 1 : print 0') == 1 ]]
    then
        echo "Add r, bioconda and conda-forge to the conda configuration."
        sudo ${CondaPATH} config --add channels r
        sudo ${CondaPATH} config --add channels bioconda
        sudo ${CondaPATH} config --add channels conda-forge
    else
        echo "Using current Conda configuration."
    fi
    {
        echo "Execute: conda create -n FASAS command ......"
        sudo ${CondaPATH} create --yes -n FASAS trimmomatic cutadapt=2.3 r-ggplot2 pigz=2.3.4 bowtie2=2.3.5 r-optparse r-pheatmap perl=5.26.2 perl-file-which perl-parallel-forkmanager
        echo "Execute: conda install --no-deps -n FASAS blast ......"
        sudo ${CondaPATH} install --yes --no-deps -n FASAS blast
    } || {
        echo "An error occurred during the execution of the conda command."
        echo "If you encounter a network connection type error, please try a few more times or modify the configuration of conda."
    }
fi

echo "Finally, install a Perl module: Switch."
unset PERL5LIB
{
    ACTIVATE_PATH=$(conda info --base)'/bin/activate'
    conda activate FASAS 2> /dev/null || source ${ACTIVATE_PATH} FASAS 2> /dev/null
    TrueInstallAdd=$(perl -e 'foreach my $every_ (@INC){ if($every_ =~ /conda/){ print $every_ and last; } }')
    NowTime=$(date +'%Y%m%d%H%M%S')
    UserTEMPFolder=${USER}'-cpanm-'${NowTime}
    cpanLog=$(cpanm Switch -l /tmp/${UserTEMPFolder})
    if [[ $(echo ${cpanLog} | perl -ne '/Successfully installed/ ? print 1 : print 0') != 1 ]]; then
        echo "Switch module is installed in other locations, exit."
        exit
    else
        pmFile=$(find /tmp/${UserTEMPFolder} -name Switch.pm)
        cp ${pmFile} ${TrueInstallAdd}
        echo "Please delete cpanm's temporary directory by yourself:"
        echo "/tmp/${UserTEMPFolder}"
    fi
} || {
    echo "Warning, the installation of the Switch module failed, please install it in the FASAS environment yourself!"
}

echo "Note One: Before using FASAS, please confirm that the conda program is in the PATH environment variable."
echo "Note Two: If the Perl module fails, note the location of the failed module."
echo "Done."
