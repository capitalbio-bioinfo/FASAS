# FASAS: Full-Length Amplicon Sequencing Analysis System
FASAS is a pipeline software for 16s full-length amplicon sequencing data analysis. The main focus of FASAS is from Rawdata to 16s rRNA contig. There are also additional functions implemented, such as: sample heatmap, krona pipline, a-diversity and b-diversity et.

## Installtion
1. Install Note: FASAS has passed the installation test in three systems: CentOS 6, CentOS 7 and Ubuntu 18.04.
2. An environment variable FASASHome is defined in FASAS. Run '. add_env.sh' or 'source add.env.sh' to write the FASASHome variable to ~/.bashrc and export it in the current SHELL.
Run:
```bash
    cd ${FASAS_HOME_FOLDER}
    . ./add_env.sh
```
3. FASAS requires some Perl modules and a conda environment, the installation script has been placed in the 'dep' directory. The installation process is very simple, just execute 'bash dep/create_conda_env.sh' and make a selection as prompted. The script will output 'done' at the end when the script runs successfully.
Run:
```bash
    cd ${FASAS_HOME_FOLDER}
    bash dep/create_conda_env.sh
```

## Quick Start
Usage:  
```
    bash run_analysis.sh -f|--config_file [ConfigFile] >& [logfile]
```
Run test:  
```
    cd [work folder]
    bash run_analysis.sh -f ./data/ConfigFile/test.ConfigFile >& test.log
```

## Parameter:
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

## Library Structure
Linker-Tag Library:  
>
>     R1:  CommonSeq---------UMISeq---------CommonSeq------------------UMISeq---------CommonSeq
>     R2:              ------UMISeq---------CommonSeq------------------UMISeq---------CommonSeq---------
>

R1 of Linker-tag Library corresponds to UMILibraryR1 of **LibraryInfo** file  
R2 of Linker-tag Library corresponds to UMILibraryR2 of **LibraryInfo** file

Read-Tag Library:  
>
>     R1:  16SSeq----------------------------------------------------------------------------------------
>                                     IN MOST CASES: NO OVERLAP
>     R2:  CommonSeq------UMISeq------CommonSeq------16SSeq----------------------------------------------
>

R1 of Read-tag Library corresponds to AssembleLibraryLeft of **LibraryInfo** file  
R2 of Read-tag Library corresponds to AssembleLibraryRight of **LibraryInfo** file

## Assemble Strategy
currently, we offer two assembly software support, [cap3](http://doua.prabi.fr/software/cap3) and [idba_ud](https://i.cs.hku.hk/~alse/hkubrg/projects/idba_ud/). this is based on the test. through testing, we believe that some genome assembly software is not suitable for Full-Length 16s data, such as: velvet, SOAP *denove*, SPAdes, ABySS. The reasons are as follows:
1. The software is cumbersome and requires too much resources
2. Poor assembly effect, not suitable for Full-Length 16s data

## Preliminary Taxonomy Algorith
**DATABASE**, The 16s database mainly includes [NCBI Blast](ftp://ftp.ncbi.nlm.nih.gov/blast/db/), [Silva](https://www.arb-silva.de/), [Greengene](http://greengenes.secondgenome.com/), [RDP](http://rdp.cme.msu.edu/) and [EzBioCloud](https://www.ezbiocloud.net).  
**Algorith**, Considering the speed and accuracy of calculations, the MegaBLAST program may be the best choice for species classification.

## Tools

Some commonly used tool-type scripts are placed in the "tools" directory. These scripts can help you analyze your data.  

## IDBA_UD

How to compile IDBA_ud supported by long sequence?  
After downloading the IDBA source code, modify line 102 of the ‘./src/sequence/short_sequence.h’ file and change the number before the semicolon at the end of the line to 152. Then compile the program normally

## License

![GPL v3 picture](https://www.gnu.org/graphics/gplv3-with-text-136x68.png)  

>    FASAS: Full-Length Amplicon Sequencing Analysis System  
>    Copyright (C) 2019 CapitalBio Corporation
>
>    This program is free software: you can redistribute it and/or modify  
>    it under the terms of the GNU General Public License as published by  
>    the Free Software Foundation, either version 3 of the License, or  
>    (at your option) any later version.
>
>    This program is distributed in the hope that it will be useful,  
>    but WITHOUT ANY WARRANTY; without even the implied warranty of  
>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  
>    GNU General Public License for more details.
>
>    You should have received a copy of the GNU General Public License  
>    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Copy Right

> CapitalBio Corporation Beijing  
> Biochip National Enginerring Research Center of Beijing
