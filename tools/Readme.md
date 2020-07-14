# Readme of Tools folder

## Create FASAS ConfigFile

- List_SequencingFile.pl    Organize sample information
- Create_ConfigFile.pl      Enter sample information and output N samples of Config file

In order to generate 16S-FASAS configuration files in batches, please complete the following preparations:

1. A 3-column table file includes sample name, sequencing name, and library type. [sample_and_sequencing.tsv]
2. Store all sequencing results in the same directory.

After determining that the preparation is complete, use the following command to batch generate the configuration file:

```shell
    perl List_SequencingFile.pl --input_file sample_and_sequencing.tsv --file_address [path of sequencing result] --split_word _ --output_file sample_information.tsv
    perl Create_ContigFile.pl --input_tsv sample_information.tsv --output_dir [your contigfile folder]
```

## Draw the cumulative curve of Contig length frequency

- Contig_Len2Num.pl    Count the number of contig lengths and plot the cumulative curve of length

## Statistical library QC information

- StatSeqNumber_number_tsv.pl    Count the QC values of the linker-tag and read-tag library sequences
- StatUMIandSeq_number_tsv.pl    Count the values of important steps in the 16S-FASAS program analysis process

## Combine multiple Blast annotation results

- Bind_TaxonomyTable.pl    Bind the Same Level Annotation Files in a 16S-FASAS Project
- Cut_LevelAnnotate.pl     Annotation Information Simplification Program

## Statistical number of taxonomy rank

- Create_TaxonomyRank.pl    Count the number of annotation at each rank of each sample

## Run ktImportText

- ktImportText.pl    ktImportText Workflow Program

## Update NCBI Blast 16S database

- Create_NCBI_BLAST_DB.pl    Get the latest 16S BLAST database from the NCBI FTP server and re-establish Taxonomy information

> *This is an unfinished job because recreating annotation information requires complex filters.

## Organize PDF files for output in the project

- Sort_out_PDF.pl    Rename and copy PDF files from a 16S-FASAS project to the same directory

## Other Tools

- Reverse_Complement.pl    DNA sequence reverse complement script
- Length_distributed.pl	   Sequence length statistics script
- Base_balance.pl          Statistical script for different base numbers
- GC_ratio.pl              A multi-threaded statistical GC ratio PERL script

