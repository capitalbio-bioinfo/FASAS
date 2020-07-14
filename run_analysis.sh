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

#check bash shell
#DBUG log:
#   zsh and ksh echo a error
#   csh report a error
QScriptEXE=$(basename `readlink /proc/$$/exe`)
( test ${QScriptEXE} = 'bash' ) || { echo -e "\033[31mError: It is recommended to use bash\033[0m"; exit 1; }

#script settings
set -o errexit

#set dafult variable
QCondaENVName='FASAS'

##########################################################################################################################
#                                                 Defind Functions                                                       #
##########################################################################################################################
function GetHelp() {
    echo "
Full-length Amplicon Sequencing Analysis System
Usage:
	bash $0 -f|--config_file [config file] > [logfile]
Run test:
        bash $0 -f ${FASASHome}/data/test.ConfigFile > test.log
Config File Format:
    '#' is commend line
    [Parameter] = [Data]
Config File Parameter:
    SampleName              strings    Sample name                               [none]
    WorkFolder              strings    Working directory                         [none]
    LibraryInfo             strings    A text file describing library design     [126_length_library_info.txt]
    LinkedLibraryR1         strings    Forward sequence of linked-tag library    [none]
    LinkedLibraryR2         strings    Reverse sequence of linked-tag library    [none]
    ReadLibraryR1           strings    Forward sequence of Read-tag library      [none]
    ReadLibraryR2           strings    Reverse sequence of Read-tag library      [none]
    CoverageDatabase        strings    Bowtie2 database of 16s rRNA              [none]
    ReferenceFasta          strings    Reference sequence, FASTA format          [NCBI 16s rRNA database]
    ReferenceTaxonomy       strings    Reference taxonomy, TSV table             [NCBI Taxonomy]
    AssembleProgram         strings    16s rRNA assembler                        [cap3 or idba_ud]
    IlluminaAdapter         strings    Sequencing adapter                        [CTGTCTCTTATACACATCT]
    ContigLength            int        Minmum contig length                      [1200]

"
    exit 0
}

function SendSignal() {
    echoerr "An error occurred during multitasking concurrency,in the $@"
    echoinfo "ready to exit."
    sleep 1
    kill $$
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

function TrapAction() {
    echoinfo "The program receives an $1 signal."
    echoinfo "The main program begins the harvesting sub-process..."
    QSubProcess=($(jobs -p))
    for PEverySubProcess in ${QSubProcess[@]}
    do
        PSubSubProcess=($(pstree -p ${PEverySubProcess} | perl -F'' -ne '@sub = $_ =~ /\((\d+)\)/g;@sub = reverse @sub;print join(" ",@sub);'))
        for PEverySubSubProcess in ${PSubSubProcess[@]}
        do
            kill -s SIGTERM ${PEverySubSubProcess} 2> /dev/null && echoinfo "SIGTERM signal has been sent to sub-process ${PEverySubSubProcess}."
        done
    done
    echoinfo "Exit."
    exit 1
}

function CPUSource() {
    local CPUSourceN=$(perl -e 'use apackages; print CheckCPU();')
    echo ${CPUSourceN}
}
##########################################################################################################################

#trap
trap 'TrapAction "Received signal SIGHUP"' SIGHUP
trap 'TrapAction "Received signal SIGINT"' SIGINT
trap 'TrapAction "Received signal SIGTERM"' SIGTERM

##if opt == 0
test $# -eq 0 && { echoerr "\033[31mThe parameter is empty, try bash $0 -h\033[0m"; exit 1; }

##get options
ConfigFile=
while true; do
    case "${1:-''}" in
        -f | --config_file) ConfigFile=$2; shift 2;;
        -h | -help |--help) GetHelp; shift 1;;
        --) shift ; break ;;
        *) break;;
    esac
done

[ ! ${ConfigFile} ] && { echoerr "\033[31mGetOptions: -f | --config_file is necessary!\033[0m"; exit 1; }

echoinfo "Active Conda environment..."
#source conda environment
#maybe, you can also use conda activate
conda activate ${QCondaENVName} 2>/dev/null || source activate ${QCondaENVName}
#check conda environment
NowCondaENV=$(conda info | grep -i 'active environment' | awk '{print $4}')
if [[ ${NowCondaENV} =~ ${QCondaENVName} ]]; then
    echoinfo "Conda ENV ${QCondaENVName} already activated!"
else
    echoerr "Conda ENV is ${NowCondaENV}, activation failed"
    echoinfo "Exit."
    exit 1
fi

#check env variable
echoinfo "Check environment variable..."
##FASASHome
if [[ ${FASASHome} =~ ^$ ]]; then
    echoerr "Environment variable FASASHome is empty."
    echoinfo "Please check the runtime environment using the env command."
    echoinfo "Exit."
    exit 1
else
    if [ ! -d ${FASASHome} ]; then
        echoerr "Environment variable FASASHome is not pointing to a folder."
        echoinfo "Check path ${FASASHome}"
        echoinfo "Exit."
        exit 1
    fi
fi
##PERL5LIB
unset PERL5LIB
export PERL5LIB=${FASASHome}'/perl/'
perl -e 'use apackages;' &> /dev/null
if [ $? != 0 ]; then
    echoerr "The value of PERL5LIB is incorrect!"
    echoinfo "Exit."
    exit 1
fi 

#analysis config file and get parameter
#SampleName WorkFolder LibraryInfo LinkedLibraryR1 LinkedLibraryR2 ReadLibraryR1 ReadLibraryR2 CoverageDatabase ReferenceFasta ReferenceTaxonomy AssembleProgram IlluminaAdapter ContigLength
{
    QParameter=($(perl ${FASASHome}/perl/00.analysis_parameter.pl --input_file ${ConfigFile}))
} || {
    if [[ ${QParameter[@]} =~ \<Error\> ]]; then
        echo ${QParameter[@]}
        echoerr "Failed to Analyze Config File, return value: $?"
        echoinfo "Exit."
        exit 1
    fi
}

#check_database
##bowtie2 database
Bowtie2IndexFileName=${QParameter[7]}'.*.bt2'
Bowtie2IndexFileNumber=$(ls ${Bowtie2IndexFileName} | wc -l)
if [ ${Bowtie2IndexFileNumber} -gt 0 ]; then
    echoinfo "Bowtie2 Database Exists"
else
    echoerr "Bowtie2 Database Error\nCheck PATH: ${QParameter[7]}"
fi

##main database
##BUG
MainIndexFolder=$(dirname ${QParameter[8]})
MainIndexFileName=${MainIndexFolder}'/*.idx'
if [ -e ${QParameter[8]} ]; then
    MainIndexFileNumber=$(ls ${MainIndexFileName} | wc -l)
    if [ ${MainIndexFileNumber} -gt 0 ]; then
        echoinfo "Main Database Index Exists"
    else
        echoerr "Main Database Index ERROR\nCheck PATH: ${QParameter[8]}"
    fi
else
    echoerr "Main Database File ERROR\nCheck PATH: ${QParameter[8]}"
fi

##main taxonomy file
test -e ${QParameter[9]} || echoerr "Main Taxonomy File ERROR\n Check PATH: ${QParameter[9]}"

#test program of assemble
which ${QParameter[10]} &> /dev/null || { echoerr "Assemble Program ${QParameter[10]} not found!"; echoinfo "Exit."; exit 1; }

echoinfo "Check input file and test size..."
#test file exists and read
test -r ${QParameter[3]} || { echoerr "${QParameter[3]} cannot be read!"; echoinfo "Exit."; exit 1; }
test -r ${QParameter[4]} || { echoerr "${QParameter[4]} cannot be read!"; echoinfo "Exit."; exit 1; }
test -r ${QParameter[5]} || { echoerr "${QParameter[5]} cannot be read!"; echoinfo "Exit."; exit 1; }
test -r ${QParameter[6]} || { echoerr "${QParameter[6]} cannot be read!"; echoinfo "Exit."; exit 1; }
test -r ${QParameter[8]} || { echoerr "${QParameter[8]} cannot be read!"; echoinfo "Exit."; exit 1; }
test -r ${QParameter[9]} || { echoerr "${QParameter[9]} cannot be read!"; echoinfo "Exit."; exit 1; }
#test contiglength
if [[ ! ${QParameter[12]} =~ ^[0-9]+$ ]]; then
    echoerr "${QParameter[12]} is not a number"
    echoinfo "Exit."
    exit 1
fi

#test file size
QLinkedSize=$(du -s ${QParameter[3]} | awk '{print $1}')
QReadSize=$(du -s ${QParameter[5]} | awk '{print $1}')
test ${QLinkedSize} -gt ${QReadSize} && { echowarn "Size of linked-tag library is bigger than Read-tag library"; }

echoinfo "Get the absolute path of the input file..."
##convert relative path to absolute path
[[ ! ${QParameter[3]} =~ ^/ ]] && QParameter[3]=$(readlink -f ${QParameter[3]})
[[ ! ${QParameter[4]} =~ ^/ ]] && QParameter[4]=$(readlink -f ${QParameter[4]})
[[ ! ${QParameter[5]} =~ ^/ ]] && QParameter[5]=$(readlink -f ${QParameter[5]})
[[ ! ${QParameter[6]} =~ ^/ ]] && QParameter[6]=$(readlink -f ${QParameter[6]})

#get work dir and cd work dir
if [[ -d ${QParameter[1]} ]]; then
    cd ${QParameter[1]}
    if [ $? != 0 ]; then
        echoerr "Work directory creation failed or no permission to access, return value $?."
        echoinfo "Exit."
        exit 1
    fi
else
    mkdir -p ${QParameter[1]} 2>/dev/null && cd ${QParameter[1]}
    if [ $? != 0 ]; then
        echoerr "Work directory creation failed or no permission to access, return value $?."
        echoinfo "Exit."
        exit 1
    fi
fi

################################################################################################################
####################                                  Main                                  ####################
####################                    step0-analysis_library_structure                    ####################
################################################################################################################

echoinfo "Begin to analysis Library..."
#QLibraryInfo Structure:
#UMILength LinkedLibraryR1P LinkedLibraryR1M LinkedLibraryR1A LinkedLibraryR2P LinkedLibraryR2M LinkedLibraryR2A ReadLibraryLeftP ReadLibraryLeftA ReadLibraryRightP ReadLibraryRightA
{
    QLibraryInfo=($(perl ${FASASHome}/perl/00.analysis_library.pl -i ${QParameter[2]}))
} || {
    echoerr "Failed to Analysis Library Info File, return value: $?."
    echoinfo "Exit."
    exit 1
}

echoinfo "Calculate library related length."
LinkedLibraryLength=$((${#QLibraryInfo[1]}+${#QLibraryInfo[2]}+${#QLibraryInfo[3]}+${QLibraryInfo[0]}*2))

echoinfo "Create soft connection for input files."
NowDirName=$(pwd)
for EveryInputFastq in ${QParameter[3]} ${QParameter[4]} ${QParameter[5]} ${QParameter[6]}; do
    InputDirName=$(dirname ${EveryInputFastq})
    InputBaseName=$(basename ${EveryInputFastq})
    if [[ (${NowDirName} != ${InputDirName}) && (! -e ${InputBaseName})]]; then
        ln -s ${EveryInputFastq} .
    fi
done

################################################################################################################
####################                          step1-split_library                           ####################
################################################################################################################
echoinfo "STEP1 split_library START!"
mkdir step1-split_library
cd step1-split_library
PUMIForwardFastqName=$(basename ${QParameter[3]})
PUMIReverseFastqName=$(basename ${QParameter[4]})
PDataForwardFastqName=$(basename ${QParameter[5]})
PDataReverseFastqName=$(basename ${QParameter[6]})
{
    perl ${FASASHome}/perl/01.split_linked_library.pl -f ../${PUMIForwardFastqName} \
        -r ../${PUMIReverseFastqName} \
        --umir1p ${QLibraryInfo[1]} -m 3 || SendSignal 'step1. script: 01.Split_UMI_Library.pl' &
    perl ${FASASHome}/perl/01.split_read_library.pl -f ../${PDataForwardFastqName} \
        -r ../${PDataReverseFastqName} \
        --alp ${QLibraryInfo[7]} --arp ${QLibraryInfo[9]} -m 3 || SendSignal 'step1. script: 01.Split_Assemble_Library.pl' &
}
wait
PStep1Work=$(pwd)
PReadLog=${PStep1Work}'/Read-tag_library.log'
PLinkedLog=${PStep1Work}'/linked-tag_library.log'
echoinfo "lead-tag library statistcs are written to ${PReadLog}"
echoinfo "linked-tag library statistcs are written to ${PLinkedLog}"

#depend log file
PMatchReadRatio=$(head -n1 Read-tag_library.log|perl -ne 'if($_ =~ /(\d+\.\d+)%\./){print $1}')
PMatchLinkedRatio=$(head -n1 linked-tag_library.log|perl -ne 'if($_ =~ /(\d+\.\d+)%\./){print $1}')
PTempNum=$(echo ${PMatchLinkedRatio} | perl -ne '$_ < 0.25 ? print 1 : print 0')
test ${PTempNum} == 1 && echowarn "Too many linked-tag library sequences do not match. Ratio: ${PMatchLinkedRatio}."
PtempNum=$(echo ${PMatchReadRatio} | perl -ne '$_ < 0.25 ? print 1 : print 0')
test ${PTempNum} == 1 && echowarn "Too many Read-tag library sequences do not match. Ratio: ${PMatchReadRatio}."

echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                     step2-trim_linked-tag_library                      ####################
################################################################################################################
echoinfo "STEP2 trim_linked-tag_library START!"
mkdir step2-trim_linked-tag_library
cd step2-trim_linked-tag_library
ln -s ../step1-split_library/linked-tag_R1.fastq.gz .
ln -s ../step1-split_library/linked-tag_R2.fastq.gz .

trimmomatic PE -threads 10 -phred33 linked-tag_R1.fastq.gz linked-tag_R2.fastq.gz \
    linked_R1_clean.fastq linked_R1_unclean.fastq.gz \
    linked_R2_clean.fastq linked_R2_unclean.fastq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:36:20 MINLEN:75 &> run_trimmomatic.log
echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                     step3-merge_linked-tag_library                     ####################
################################################################################################################
echoinfo "STEP3 merge_linked-tag_library START!"
mkdir step3-merge_linked-tag_library
cd step3-merge_linked-tag_library
ln -s ../step2-trim_linked-tag_library/linked_R1_clean.fastq .
ln -s ../step2-trim_linked-tag_library/linked_R2_clean.fastq .
#XORRO not suppect degenerate base.
#XORRO not suppect gzip format.
#XORRO only output the sequence of the overlap, which meets the requirements of the linked-tag library.
perl ${FASASHome}/bin/xorro_wrapper.pl -i1 linked_R1_clean.fastq -i2 linked_R2_clean.fastq \
    -o linked_merge.assembled.fastq -k T

#length of good linked-tag library sequence < ( library length + 24 )
PLinkedCutoff=$((${LinkedLibraryLength}+24))
#stat sequence length of the linked-tag library
cat linked_merge.assembled.fastq | perl -ne '$seq=<>;$mid=<>;$qual=<>;chomp $seq;print length $seq,"\n";' > linked_sequence_length.txt
Rscript ${FASASHome}/r/step3-plotting_linked-tag_length.r linked_sequence_length.txt ${PLinkedCutoff} linked-tag_length_frequency.pdf 2> /dev/null
#delete some sequence that sequence length gt ${PLinkedCutoff}
perl ${FASASHome}/perl/03.filter_linked-tag_sequence.pl --upper_limit ${PLinkedCutoff} --input_file linked_merge.assembled.fastq --output_file linked_merge.assembled.fasta

#correct UMI library
#require:
#similarity > 80%, exist umi and umi length le fixed length
#This a Bad Alignmnet Result, Example:
#TGAGCCAKGATCAAACTCTTAGATCGCNNNNNNNNNNNNNNATGGATGAGTCTGGGTGGAGAGGGGGGCAAAGATGAAGATNNNNNNNNNNNNNNCGTACTAGTACGGYTACCTTGTTACGACTT
#TGAGCCAGGATCAAACTCTTAGATCGCCTACGAGAGGGCACATGGATGAGTCTGGGTGGAG-------------T-A------------------C-TA---GTACGGTTACCTTGTTACGACTT
perl ${FASASHome}/perl/03.correct_linked_sequence.pl --input_fasta linked_merge.assembled.fasta \
    --output_file linked_correct.fasta --ur1p ${QLibraryInfo[1]} --ur1m ${QLibraryInfo[2]} \
    --ur1a ${QLibraryInfo[3]} --umi_length ${QLibraryInfo[0]} \
    --umi_toolong_log linked_too_long.fasta --low_similarity_log linked_correct.log
echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                    step4-stat_linked-tag_library                       ####################
################################################################################################################
echoinfo "STEP4 stat_linked_tag_library START!"
mkdir step4-stat_linked-tag_library
cd step4-stat_linked-tag_library
ln -s ../step3-merge_linked-tag_library/linked_correct.fasta .
perl ${FASASHome}/perl/04.create_UMI_table.pl -f linked_correct.fasta \
    -l left_UMItable.tsv -r right_UMItable.tsv \
    --ur1p ${QLibraryInfo[1]} --ur1m ${QLibraryInfo[2]} \
    --ur1a ${QLibraryInfo[3]} --umi_length ${QLibraryInfo[0]} --log_file create_UMI_table.log

{
    #data pre-thread
    perl ${FASASHome}/perl/04.output_for_plot.pl -i left_UMItable.tsv -o left_UMI_frequency.tsv || SendSignal 'step4. script: 04.OutputForPlot.pl LEFT' &
    perl ${FASASHome}/perl/04.output_for_plot.pl -i right_UMItable.tsv -o right_UMI_frequency.tsv || SendSignal 'step4. script: 04.OutputForPlot.pl RIGHT' &
}
wait
Rscript ${FASASHome}/r/step4-plotting_UMI_frequency.r left_UMI_frequency.tsv right_UMI_frequency.tsv UMI_frequency.pdf 2> /dev/null
echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                      step5-trim_read-tag_library                       ####################
################################################################################################################
echoinfo "STEP5 trim_read-tag_library START!"
mkdir step5-trim_read-tag_library
cd step5-trim_read-tag_library
ln -s ../step1-split_library/read-tag_left_R1.fastq.gz .
ln -s ../step1-split_library/read-tag_left_R2.fastq.gz .
ln -s ../step1-split_library/read-tag_right_R1.fastq.gz .
ln -s ../step1-split_library/read-tag_right_R1.fastq.gz .
#max_threads option is 10, actually, 6 is top
{
    perl ${FASASHome}/perl/05.takeout_UMI.pl -f read-tag_left_R1.fastq.gz -r read-tag_left_R2.fastq.gz \
        -t l --aleftp ${QLibraryInfo[7]} --alefta ${QLibraryInfo[8]} \
        --arightp ${QLibraryInfo[9]} --arighta ${QLibraryInfo[10]} \
        --umi_length ${QLibraryInfo[0]} --max_threads 10 || SendSignal 'step5. script: 05.TakeOut_UMI.pl LEFT' &
    perl ${FASASHome}/perl/05.takeout_UMI.pl -f read-tag_right_R1.fastq.gz -r read-tag_right_R1.fastq.gz \
        -t r --aleftp ${QLibraryInfo[7]} --alefta ${QLibraryInfo[8]} \
        --arightp ${QLibraryInfo[9]} --arighta ${QLibraryInfo[10]} \
        --umi_length ${QLibraryInfo[0]} --max_threads 10 || SendSignal 'step5. script: 05.TakeOut_UMI.pl RIGHT' &
}
wait

echoinfo "Run cutadapt..."
AvialableCPU=$(CPUSource)
test ${AvialableCPU} > 20 && AvialableCPU=20
##################################################################################
#cutadapt=2.3 python api in conda
#${QLibraryInfo[1]} is a 16s primer, at the front of 16srRNA
#${QLibraryInfo[4]} is a 16s primer, at the end of 16srRNA
#-O overlap >= 6
#-m min_length >= 50
#-q quality_cutoff PE >= 20
#-j threads == 20
##################################################################################

cutadapt -a ${QLibraryInfo[1]} -a ${QLibraryInfo[7]} -g ${QLibraryInfo[10]} \
    -A ${QLibraryInfo[4]} -G ${QParameter[11]} \
    -o read-tag_left_R1_two.fastq.gz -p read-tag_left_R2_two.fastq.gz read-tag_left_R1_one.fastq.gz read-tag_left_R2_one.fastq.gz \
    -O 6 -m 50 -q 20,20 \
    -j ${AvialableCPU} 1> cutadapt_left.log 2> /dev/null || SendSignal 'step5. script: cutadapt LEFT'
cutadapt -a ${QLibraryInfo[4]} -a ${QLibraryInfo[9]} -g ${QLibraryInfo[8]} \
    -A ${QLibraryInfo[1]} -G ${QParameter[11]} \
    -o read-tag_right_R1_two.fastq.gz -p read-tag_right_R2_two.fastq.gz read-tag_right_R1_one.fastq.gz read-tag_right_R2_one.fastq.gz \
    -O 6 -m 50 -q 20,20 \
    -j ${AvialableCPU} 1> cutadapt_right.log 2> /dev/null || SendSignal 'step5. script: cutadapt RIGHT'

echoinfo "Filter Sequence..."
{
    #Goal: Delete ploy END.
    #--------------------------LEFT 16s rRNA-------------------------
    #                                   <-- R1 (sequencing direction)
    perl ${FASASHome}/perl/05.filter_R1.pl --input_r1 read-tag_left_R1_two.fastq.gz \
        --input_r2 read-tag_left_R2_two.fastq.gz \
        --output_r1 read-tag_left_R1_three.fastq.gz --output_r2 read-tag_left_R2_three.fastq.gz \
        --format fq --type LEFT || SendSignal 'step5. script: 05.Filter_R1.pl LEFT' &
    perl ${FASASHome}/perl/05.filter_R1.pl --input_r1 read-tag_right_R1_two.fastq.gz \
        --input_r2 read-tag_right_R2_two.fastq.gz \
        --output_r1 read-tag_right_R1_three.fastq.gz --output_r2 read-tag_right_R2_three.fastq.gz \
        --format fq --type RIGHT || SendSignal 'step5. script: 05.Filter_R1.pl RIGHT' &
}
wait
echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                      step6-stat_read-tag_library                       ####################
################################################################################################################
echoinfo "STEP6 stat_read-tag_library START!"
mkdir step6-stat_read-tag_library
cd step6-stat_read-tag_library
ln -s ../step5-trim_link_library/read-tag_left_R2_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_right_R2_three.fastq.gz .
{
	perl ${FASASHome}/perl/06.stat_read_library_reads.pl -r read-tag_left_R2_three.fastq.gz -t l \
        --output_number left_UMI2number.tsv \
        --output_seqname left_UMI2seqname.tsv || SendSignal 'step6. script: 06.StatAssembleRead.pl LEFT' &
	perl ${FASASHome}/perl/06.stat_read_library_reads.pl -r read-tag_right_R2_three.fastq.gz -t r \
        --output_number right_UMI2number.tsv \
        --output_seqname right_UMI2seqname.tsv || SendSignal 'step6. script: 06.StatAssembleRead.pl RIGHT' &
}
wait
echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                         step7-determine_sequence                       ####################
################################################################################################################
echoinfo "STEP7 determine_sequence START!"
mkdir step7-determine_sequence
cd step7-determine_sequence
ln -s ../step4-stat_UMI_library/left_UMItable.tsv .
ln -s ../step4-stat_UMI_library/right_UMItable.tsv .
ln -s ../step5-trim_link_library/read-tag_left_R1_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_left_R2_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_right_R1_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_right_R2_three.fastq.gz .
ln -s ../step6-stat_link_library/left_UMI2seqname.tsv .
ln -s ../step6-stat_link_library/right_UMI2seqname.tsv .
{
    perl ${FASASHome}/perl/07.16s_distributor.pl --left_input_umi2name left_UMI2seqname.tsv \
        --right_input_umi2name right_UMI2seqname.tsv \
        --left_input_umitable left_UMItable.tsv --right_input_umitable right_UMItable.tsv \
        --left_input_r1fastq read-tag_left_R1_three.fastq.gz --right_input_r1fastq read-tag_right_R1_three.fastq.gz \
        --left_input_r2fastq read-tag_left_R2_three.fastq.gz --right_input_r2fastq read-tag_right_R2_three.fastq.gz --filter_read_number 25
} || {
    echoerr "Failed to Split Assemble Library, return value: $?."
    echoinfo "Exit."
    exit 1
}

PAllSeqSize=$(du -s All_seq_r1.fastq.gz|awk '{print $1}')
if [ ! ${PAllSeqSize} > 40 ]; then
    echoerr "File All_seq_r1.fastq has a size of zero!"
fi

#check_seq_cover
cd ../
mkdir check_read-tag_coverage
cd check_read-tag_coverage
ln -s ../step5-trim_link_library/read-tag_left_R1_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_left_R2_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_right_R1_three.fastq.gz .
ln -s ../step5-trim_link_library/read-tag_right_R2_three.fastq.gz .
#sampling
{
    perl ${FASASHome}/perl/07.sampling_reads.pl --input_r1 read-tag_left_R1_three.fastq.gz \
        --input_r2 read-tag_left_R2_three.fastq.gz \
        --output_r1 sampling_left_100000_R1.fasta.gz --output_r2 sampling_left_100000_R2.fasta.gz \
        --sampling_number 100000 || SendSignal 'step7. script: 07.PE_Read_Sampling.pl LEFT' &
    perl ${FASASHome}/perl/07.sampling_reads.pl --input_r1 read-tag_right_R1_three.fastq.gz \
        --input_r2 read-tag_right_R2_three.fastq.gz \
        --output_r1 sampling_right_100000_R1.fasta.gz --output_r2 sampling_right_100000_R2.fasta.gz \
        --sampling_number 100000 || SendSignal 'step7. script: 07.PE_Read_Sampling.pl RIGHT' &
}
wait
{
    bowtie2 --sensitive-local -f sampling_left_100000_R1.fasta.gz -x ${QParameter[7]} --no-hd --no-unal  --no-sq --quiet -p 5 > sampling_left_r1.sam || SendSignal 'step7. bowtie2 left_r1' &
    bowtie2 --sensitive-local -f sampling_left_100000_R2.fasta.gz -x ${QParameter[7]} --no-hd --no-unal  --no-sq --quiet -p 5 > sampling_left_r2.sam || SendSignal 'step7. bowtie2 left_r2' &
    bowtie2 --sensitive-local -f sampling_right_100000_R1.fasta.gz -x ${QParameter[7]} --no-hd --no-unal  --no-sq --quiet -p 5 > sampling_right_r1.sam || SendSignal 'step7. bowtie2 right_r1' &
    bowtie2 --sensitive-local -f sampling_right_100000_R2.fasta.gz -x ${QParameter[7]} --no-hd --no-unal  --no-sq --quiet -p 5 > sampling_right_r2.sam || SendSignal 'step7. bowtie2 right_r2' &
}
wait
{
    perl ${FASASHome}/perl/07.read_sam.pl -input_path ./ -type r1 -o sampling_coverage_r1.cov || SendSignal 'step7. script: 07.Read_Sam.pl R1' &
    perl ${FASASHome}/perl/07.read_sam.pl -input_path ./ -type r2 -o sampling_coverage_r2.cov || SendSignal 'step7. script: 07.Read_Sam.pl R2' &
}
wait
Rscript ${FASASHome}/r/step7-plotting_reads_coverage.r sampling_coverage_r1.cov sampling_coverage_r2.cov 100k_pe_sampling_coverage.pdf 2> /dev/null
echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                        step8-parallel_assembly                         ####################
################################################################################################################
echoinfo "STEP8 parallel_assembly START!"
mkdir step8-parallel_assembly
cd step8-parallel_assembly
ln -s ../step7-split_link_library/all_sequence_R1.fastq.gz .
ln -s ../step7-split_link_library/all_sequence_R2.fastq.gz .

#BUG:
#threads options: -m ???
{
    perl ${FASASHome}/perl/08.parallel_assemble.pl -f all_sequence_R1.fastq.gz -r all_sequence_R2.fastq.gz -o contigs.fa -m 40 -a ${QParameter[10]} 2> /dev/null
} || {
    echoerr "Failed to Parallel Assembly, return value: $?."
    echoinfo "Exit."
    exit 1
}
if [ ${QParameter[10]} == 'idba_ud' ]; then
    perl ${FASASHome}/perl/08.close_gap.pl --folder_log ./idba_ud_address.txt \
        --contigs ./contigs.fa --db ${QParameter[7]} \
        --contig_limit 1200 --max_threads 20 \
        --output second_contigs.fa --log close_gap.log
fi

perl ${FASASHome}/perl/08.get_contig_length.pl --min_value 300 --max_value 1800 --input_file second_contigs.fa --output_file contig_length.log
Rscript ${FASASHome}/r/step8-plotting_contig_frequency.r contig_length.log contig_frequency.pdf 2> /dev/null

PContigsNumber=$(wc -l fix_contigs.fa|awk '{print $1}')
test ${PContigsNumber} -le 6000 && echowarn "STEP8: The number of Contig is less than 6000!"

echoinfo "=================================================="
cd - > /dev/null

################################################################################################################
####################                      step9_megablast_annotation                        ####################
################################################################################################################
echoinfo "STEP9 megablast_annotation START!"
mkdir step9_megablast_annotation
cd step9_megablast_annotation
if [ ${QParameter[10]} == 'idba_ud' ]; then
    QStep9Input='second_contigs.fa'
    ln -s ../step8-parallel_assembly/second_contigs.fa .
else
    QStep9Input='contigs.fa'
    ln -s ../step8-parallel_assembly/contigs.fa .
fi
#${QParameter[8]} need fasta format, because the usearch11 program needs
perl ${FASASHome}/perl/09.megablast_annotation.pl --work_dir ./ \
    --input_fasta ${QStep9Input} --ref_fa ${QParameter[8]} \
    --ref_tax ${QParameter[9]} --blast_db ${QParameter[8]} \
    --contig_length ${QParameter[12]} --blast_evalue 1e-20
perl ${FASASHome}/perl/09.stat_taxonmy.pl --input_file ./primer_taxonomy.tsv  --output_dir ./every_rank --sample_name ${QParameter[0]}
echoinfo "=================================================="
echoinfo "The program is successfully completed, END."
################################################################################################################
####################                             sample_report                              ####################
################################################################################################################
##NOTE:In the future, I have the opportunity, use LeTex to make a simple report.

