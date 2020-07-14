#!/usr/bin/perl
=head
    FileName: stat_reads_procedure.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    stat_reads_procedure.pl
    Copyright (C) 2019 CapitalBio Corporation

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
=cut

=function
    Input: FASAS project folder
    Output: Linked-tag library and Read-tag library statistcs result
=cut

use strict;
use warnings;
use apackages;
use Cwd qw(abs_path);
use IO::File;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my ($search_dir,$number_support,$help);
my $linked_tag = '';
my $read_tag = '';
$number_support = 25;

my $options_number = scalar @ARGV;
GetOptions (
    "search_dir=s" => \$search_dir,
    "linked_tag=s" => \$linked_tag,
    "read_tag=s" => \$read_tag,
    "number_support=i" => \$number_support,
    "help" => \$help,
) or die "Failed to Get Options\n";

if ($help or $options_number < 2){
    print "WARN: option number lt 2!\n" if $options_number < 2;
    print "Sequence Number Statistics Program
Usage:
    perl $0 --search_dir [16S-FASAS project folder] --linked_tag [umi file] --read_tag [assemble file]
Options:
    search_dir       path      
    linked_tag       strings   
    read_tag         strings   
    number_support   int    
    help             switch    
";
    exit;
}

my $UMIFile = 'step1-split_library/Linked-tag_library.log';
my $UMICleanFile = 'step3-merge_UMI_library/UMI_Correct.fasta';
my $AssembleFile = 'step1-split_library/Read-tag_library.log';
my $AssembleLeftCleanFile = 'step5-trim_link_library/LEFT_R1_StepThree.fastq.gz';
my $AssembleRightCleanFile = 'step5-trim_link_library/RIGHT_R1_StepThree.fastq.gz';
my $AssembleLeftStatFile = 'step6-stat_link_library/LEFT_UMI2Number.tsv';
my $AssembleRightStatFile = 'step6-stat_link_library/RIGHT_UMI2Number.tsv';
my $AssembleLastFile = 'step7-split_link_library/All_seq_r1.fastq.gz';

my @folder = @{GetFolder($search_dir)};


my %result;
my @sample_name;
my $sample_number = 0;

foreach my $every_sample (sort { $a cmp $b } @folder){
    #samplename
    my $sample_name = basename($every_sample);
    #get file address
    my $umi_stat_file = $every_sample.'/'.$UMIFile;
    my $assemble_stat_file = $every_sample.'/'.$AssembleFile;
    #read file
    my $umi_stat_info = ReadLogFile($umi_stat_file);
    my $assemble_stat_info = ReadLogFile($assemble_stat_file);
    ##begin stat
    #UMI
    #RawReads
    $result{$sample_name}{'RawReads'} = $1 if $umi_stat_info =~ /Total\s+Seq:\s+(\d+)/;
    #CleanReads
    my $umi_clean_file = $every_sample.'/'.$UMICleanFile;
    $result{$sample_name}{'CleanReads'} = GetFileLine($umi_clean_file) / 2;
    $result{$sample_name}{'Ratio'} = $result{$sample_name}{'CleanReads'} / $result{$sample_name}{'RawReads'} * 100;
    #Assemble
    $result{$sample_name}{'Read-tag RawReads'} = $1 if $assemble_stat_info =~ /Total\s+Seq:\s+(\d+)/;
    $result{$sample_name}{'Left RawReads'} = $1 if $assemble_stat_info =~ /LEFT:\s+(\d+)/;
    $result{$sample_name}{'Right RawReads'} = $1 if $assemble_stat_info =~ /RIGHT:\s+(\d+)/;
    $result{$sample_name}{'Left Prop.'} = $1 if $assemble_stat_info =~ /LEFT:\s+\d+,pr:\s+(\d+\.\d+)/;
    $result{$sample_name}{'Right Prop.'} = $1 if $assemble_stat_info =~ /RIGHT:\s+\d+,pr:\s+(\d+\.\d+)/;
    #clean stat
    my $assemble_left_clean_file = $every_sample.'/'.$AssembleLeftCleanFile;
    my $assemble_right_clean_file = $every_sample.'/'.$AssembleRightCleanFile;
    $result{$sample_name}{'Left CleanReads'} = GetFileLine($assemble_left_clean_file) / 4;
    $result{$sample_name}{'Right CleanReads'} = GetFileLine($assemble_right_clean_file) / 4;
    $result{$sample_name}{'Read-tag CleanReads'} = $result{$sample_name}{'Left CleanReads'} + $result{$sample_name}{'Right CleanReads'};
    $result{$sample_name}{'LeftCleanRatio'} = $result{$sample_name}{'Left CleanReads'} / $result{$sample_name}{'Left RawReads'} * 100;
    $result{$sample_name}{'RightCleanRatio'} = $result{$sample_name}{'Right CleanReads'} / $result{$sample_name}{'Right RawReads'} * 100;
    $result{$sample_name}{'Read-tagCleanRatio'} = $result{$sample_name}{'Read-tag CleanReads'} / $result{$sample_name}{'Read-tag RawReads'} * 100;
    #last stat
    my $assemble_left_seq2num_file = $every_sample.'/'.$AssembleLeftStatFile;
    my $assemble_right_seq2num_file = $every_sample.'/'.$AssembleRightStatFile;
    my $assemble_last_file = $every_sample.'/'.$AssembleLastFile;
    $result{$sample_name}{'LeftGoodNumber'} = GetGoodNumber($assemble_left_seq2num_file);
    $result{$sample_name}{'RightGoodNumber'} = GetGoodNumber($assemble_right_seq2num_file);
    $result{$sample_name}{'Last Reads'} = GetFileLine($assemble_last_file) / 4;
    $result{$sample_name}{'Clean/Last Ratio'} = $result{$sample_name}{'Last Reads'} / ( $result{$sample_name}{'Left CleanReads'} + $result{$sample_name}{'Right CleanReads'} ) * 100;
    $result{$sample_name}{'Good/Last Ratio'} = $result{$sample_name}{'Last Reads'} / ( $result{$sample_name}{'LeftGoodNumber'} + $result{$sample_name}{'RightGoodNumber'} ) * 100;

    $sample_number ++;
    push @sample_name,$sample_name
}

##Total
foreach my $every_sample (keys %result){
    foreach my $every_attr (keys %{$result{$every_sample}}){
        $result{'Total'}{$every_attr} += $result{$every_sample}{$every_attr};
    }
}
$result{'Total'}{'Read-tagCleanRatio'} = $result{'Total'}{'Read-tagCleanRatio'} / $sample_number;
$result{'Total'}{'Ratio'} = $result{'Total'}{'Ratio'} / $sample_number;
$result{'Total'}{'LeftCleanRatio'} = $result{'Total'}{'LeftCleanRatio'} / $sample_number;
$result{'Total'}{'RightCleanRatio'} = $result{'Total'}{'RightCleanRatio'} / $sample_number;
$result{'Total'}{'Left Prop.'} = $result{'Total'}{'Left Prop.'} / $sample_number;
$result{'Total'}{'Right Prop.'} = $result{'Total'}{'Right Prop.'} / $sample_number;
$result{'Total'}{'Clean/Last Ratio'} = $result{'Total'}{'Clean/Last Ratio'} / $sample_number;
$result{'Total'}{'Good/Last Ratio'} = $result{'Total'}{'Good/Last Ratio'} / $sample_number;
push @sample_name,'Total';

##format say
my @umi_array = qw(SampleID RawReads CleanReads Ratio);
my @assemble_array = ('SampleID','Left RawReads','Left Prop.','Left CleanReads','LeftCleanRatio','Right RawReads','Right Prop.','Right CleanReads','RightCleanRatio','Read-tag RawReads','Read-tag CleanReads','Read-tagCleanRatio','LeftGoodNumber','RightGoodNumber','Last Reads','Clean/Last Ratio','Good/Last Ratio');

if ($linked_tag eq ''){
    print join("\t",@umi_array),"\n";
    foreach my $every_sample (@sample_name){
        print "$every_sample";
        for (my $i = 1;$i <= $#umi_array; $i ++){
            if ($umi_array[$i] =~ /ratio|prop/i){
                printf "\t%.2f%%",$result{$every_sample}{$umi_array[$i]};
            }else{
                print "\t$result{$every_sample}{$umi_array[$i]}";
            }
        }
        print "\n";
    }
}else{
    my $umi_output_handle;
    $umi_output_handle = new IO::File ">$linked_tag" or die "$!";

    print $umi_output_handle join("\t",@umi_array),"\n";
    foreach my $every_sample (@sample_name){
        print $umi_output_handle "$every_sample";
        for (my $i = 1;$i <= $#umi_array; $i ++){
            if ($umi_array[$i] =~ /ratio|prop/i){
                printf $umi_output_handle "\t%.2f%%",$result{$every_sample}{$umi_array[$i]};
            }else{
                print $umi_output_handle "\t$result{$every_sample}{$umi_array[$i]}";
            }
        }
        print $umi_output_handle "\n";
    }
    $umi_output_handle -> close;
}

if ($read_tag eq ''){
    print join("\t",@assemble_array),"\n";
    foreach my $every_sample (@sample_name){
        print "$every_sample";
        for (my $i = 1;$i <= $#assemble_array; $i ++){
            if ($assemble_array[$i] =~ /ratio|prop/i){
                printf "\t%.2f%%",$result{$every_sample}{$assemble_array[$i]};
            }else{
                print "\t$result{$every_sample}{$assemble_array[$i]}";
            }
        }
        print "\n";
    }
}else{
    my $assemble_output_handle;
    $assemble_output_handle = new IO::File ">$read_tag" or die "$!";

    print $assemble_output_handle join("\t",@assemble_array),"\n";
    foreach my $every_sample (@sample_name){
        print $assemble_output_handle "$every_sample";
        for (my $i = 1;$i <= $#assemble_array; $i ++){
            if ($assemble_array[$i] =~ /ratio|prop/i){
                printf $assemble_output_handle "\t%.2f%%",$result{$every_sample}{$assemble_array[$i]};
            }else{
                print $assemble_output_handle "\t$result{$every_sample}{$assemble_array[$i]}";
            }
        }
        print $assemble_output_handle "\n";
    }
    $assemble_output_handle -> close;
}

print Dumper %result;

sub ReadLogFile {
    my $file = shift;
    my $file_handle = new IO::File "$file" or die "$!";
    local $/ = undef;
    my $file_info = <$file_handle>;
    return $file_info;
}

sub GetFolder {
    my $search_target = shift;
    my @search_result;
    opendir SEARCH,"$search_target" or die "$!";
    while (my $next_source = readdir SEARCH){
        next if $next_source eq '..';
        my $source_path = $search_target.'/'.$next_source;
        $source_path = abs_path($source_path);
        if (-d $source_path){
            opendir NEXTSOURCE,"$source_path" or die "$!";
            my $step_check_status = 0;
            while (my $check_source = readdir NEXTSOURCE){
                if ($check_source =~ /step/){
                    $step_check_status = 1;
                    last;
                }
            }
            closedir NEXTSOURCE;
            if ($step_check_status == 1){
                push @search_result,$source_path;
            }else{
                next;
            }
        }else{
            next;
        }
    }
    closedir SEARCH;
    return \@search_result;
}

sub GetFileLine {
    my $input_file = shift;
    my $file_line_number = 0;
    my $input_strings = gzip_support("$input_file","input",2);
    print "$input_strings\n";
    my $input_file_handle = new IO::File "$input_strings" or die "$input_file\n$!";
    while (<$input_file_handle>){
        $file_line_number ++;
    }
    $input_file_handle -> close;
    return $file_line_number;
}

sub GetGoodNumber {
    my $file = shift;
    my $input_handle = new IO::File "$file" or echo_error("Stat","$!");
    my $good_number = 0;
    while (<$input_handle>){
        chomp;
        my @line = split "\t",$_;
        $good_number += $line[1] if $line[1] >= 25;
    }
    $input_handle -> close;
    return $good_number;
}
