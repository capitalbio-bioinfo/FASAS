# 16S-FASAS: 16S rRNA full-Length Amplicon Sequencing Analysis System

## Introduction
16S-FASAS is a full-length 16S amplicon sequencing data analysis system which contains collections of modules such as data quality control, sequence demultiplexing, parallel assembly, taxonomy annotation and so on. 
![Overall workflow of 16S-FASAS](https://github.com/capitalbio-bioinfo/FASAS/blob/master/data/project_png/figure1.jpg)

## Installtion
1. Clone the repository to your local path. 
Run:
```bash
    mkdir /PATH/TO/FASASHome/
    git clone https://github.com/capitalbio-bioinfo/FASAS.git /PATH/TO/FASASHome/
```
2. An environment variable FASASHome is defined in 16S-FASAS. Run 'bash add_env.sh' or 'source add.env.sh' to write the FASASHome variable to ~/.bashrc and export it in the current SHELL.
Run:
```bash
    cd /PATH/TO/FASASHome/
    bash ./add_env.sh
```

3. 16S-FASAS requires some Perl modules and a conda environment, the installation script has been placed in the 'dep' directory. The installation process is very simple, just execute 'bash dep/create_conda_env.sh' and make a selection as prompted. The script will output 'done' at the end when the script runs successfully.
Run:
```bash
   bash dep/create_conda_env.sh
```

4. To assemble paired-end reads with longer read length, you have to re-compile IDBA_ud from the source as follows. Modify the constant "kMaxShortSequence" from 128 to 157 (or more) in /src/sequence/short_sequence.h:

**before:**
 ```
    102     static const uint32_t kMaxShortSequence = 128;  
 ```
**after:**
```
    102     static const uint32_t kMaxShortSequence = 157;
```
```bash
    cd /PATH/TO/idba/idba-master
    configure & make!
```
 
## Usage
All parameters are specified in the analyzer configuration file.
```bash
    bash run_analysis.sh -f|--config_file [ConfigFile] >& [logfile]
```

## Configuration parameters
```bash
SampleName              strings    Sample name                               [none]
WorkFolder              strings    Working directory                         [none]
LibraryInfo             strings    A text file describing library design     [126_length_library_info.txt]
LinkedLibraryR1         strings    Forward sequence of linked-tag library    [none]
LinkedLibraryR2         strings    Reverse sequence of linked-tag library    [none]
ReadLibraryR1           strings    Forward sequence of Read-tag library      [none]
ReadLibraryR2           strings    Reverse sequence of Read-tag library      [none]
CoverageDatabase        strings    Bowtie2 database of 16s rRNA              [none]
ReferenceFasta          strings    Reference sequence, FASTA format          [none]
ReferenceTaxonomy       strings    Reference taxonomy, TSV table             [none]
AssembleProgram         strings    16s rRNA assembler                        [cap3 or idba_ud]
IlluminaAdapter         strings    Sequencing adapter                        [none]
ContigLength            int        Minmum contig length                      [1200]
```

## Test
Run the following codes for test:
```bash
    #build test database
    bash ${FASASHome}/data/example_data/database/build_database.sh ${FASASHome}/data/example_data/database/mini_fulllength.fasta
    #create config file
    cd [work folder]
    bash ${FASASHome}/data/example_data/create_test_config.sh > ./test.configfile
    #run test analysis
    bash ${FASASHome}/run_analysis.sh -f ./test.configfile >& test.log
```


## Database
 Database could be built based on the [NCBI Taxonomy database](https://ftp.ncbi.nih.gov/pub/taxonomy/), [Silva](https://www.arb-silva.de/), [Greengene](http://greengenes.secondgenome.com/), [RDP](http://rdp.cme.msu.edu/) and [EzBioCloud](https://www.ezbiocloud.net).  

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
