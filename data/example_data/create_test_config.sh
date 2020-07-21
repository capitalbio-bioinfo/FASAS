#!/bin/bash
#Create Test Config File
#Usage:
#bash $0 > test.configfile

CONFIG_DATA="cat"
"${CONFIG_DATA}" <<EOF
SampleName = test_sample
WorkFolder = test_sample
LibraryInfo = ${FASASHome}/data/library_info/126_length_library_info.txt
LinkedLibraryR1 = ${FASASHome}/data/example_data/raw_data/20200509_16S_63_R1.fq.gz
LinkedLibraryR2 = ${FASASHome}/data/example_data/raw_data/20200509_16S_63_R2.fq.gz
ReadLibraryR1 = ${FASASHome}/data/example_data/raw_data/20200509_16S_26_R1.fq.gz
ReadLibraryR2 = ${FASASHome}/data/example_data/raw_data/20200509_16S_26_R2.fq.gz
CoverageDatabase = ${FASASHome}/data/example_data/database/mini_fulllength.fasta
ReferenceFasta = ${FASASHome}/data/example_data/database/mini_fulllength.fasta
ReferenceTaxonomy = ${FASASHome}/data/example_data/database/mini_taxonomy.txt
AssembleProgram = idba_ud
IlluminaAdapter = CTGTCTCTTATACACATCT
ContigLength = 1200
EOF
