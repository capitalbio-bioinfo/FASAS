#!/usr/bin/perl
=head
    FileName: stat_FASAS_procedure.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    stat_FASAS_procedure.pl
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
    Output: FASAS procedure statistcs result
=cut

use strict;
use warnings;
use apackages;
use Cwd qw(abs_path);
use IO::File;
use List::Util qw(sum);
use Getopt::Long;
use File::Basename;

my ($search_dir,$output_file,$long_contig_length,$reads_support_number,$help);

#default options
$long_contig_length = 1200;

my $options_number = scalar @ARGV;
GetOptions (
    "search_dir=s" => \$search_dir,
    "output_file=s" => \$output_file,
    "long_contig_length=i" => \$long_contig_length,
    "reads_support_number=i" => \$reads_support_number,
    "help" => \$help,
) or die "Failed to Get Options\n";

if ($help or $options_number < 4){
    print "WARN: option number lt 4!\n\n" if $options_number < 4;
    print "Count Program for FASAS Analysis Procedure
Usage:
    perl $0 -s [search_folder] -o [output_tsv_file]
Options:
    search_dir            path       Input folder in which to search the project directory [required]
    output_file           strings    Output file, tsv format [required]
    long_contig_length    int        As a threshold for effective Contig [1200]
    reads_support_number  int        Lowest sequence support for each UMI in the read-tag library [require]
    help                  switch     Get this message [off]

";
    exit;
}

if (defined $reads_support_number){
    if ($reads_support_number <= 10){
        echo_warn("Stat_FASAS","Support number of reads is $reads_support_number, do we really use this threshold?");
    }
}else{
    echo_error("Stat_FASAS","Need reads_support_number options!");
}

#search folder
my @folder = @{GetFolder($search_dir)};
die "StatUMIandSeq_number_tsv: No project directory found in folder $search_dir !\n" if scalar @folder == 0;
#defined file name
my $LtoRFile = 'step4-stat_linked-tag_library/StatUMI_LtoR.tsv';
my $RtoLFile = 'step4-stat_linked-tag_library/StatUMI_RtoL.tsv';
my $LEFTSeqFile = 'step6-stat_read-tag_library/LEFT_UMI2Number.tsv';
my $RIGHTSeqFile = 'step6-stat_read-tag_library/RIGHT_UMI2Number.tsv';
my $Split16SFile = 'step7-split_link_library/All_seq_r1.fastq.gz';
my $Contig = 'step8-parallel_assembly/fix_contigs.fa';
my $ClearContig = 'step9-megablast_annotation/contigs_clear.fa';

#get poutput handle
my $output_handle;
$output_handle = new IO::File ">$output_file" or die "$!";

#header
print $output_handle "SampleName\tLEFTUMIN.\tRIGHTUMIN.\tLEFTGood16SN.\tLEFTGoodSeqN.\tLEFTBad16SN.\tLEFTBadSeqN.\tRIGHTGood16SN.\tRIGHTGoodSeqN.\tRIGHTBad16SN.\tRIGHTBadSeqN.\t16SNumber\tRawContigNumber\tContigN50\tRawContigNumber(length >= $long_contig_length)\tClearContigNumber\n";
foreach my $every_dir (sort { $a cmp $b } @folder){
    #dirname
    my $dirname = basename($every_dir);
    #linked-tag library
    my $LEFT_umi_file = $every_dir.'/'.$LtoRFile;
    my $UMI_LEFT_number;
    if (-s $LEFT_umi_file){
        $UMI_LEFT_number = GetFileLine($LEFT_umi_file);
    }else{
        $UMI_LEFT_number = 0;
    }
    my $RIGHT_umi_file = $every_dir.'/'.$RtoLFile;
    my $UMI_RIGHT_number;
    if (-s $RIGHT_umi_file){
        $UMI_RIGHT_number = GetFileLine($RIGHT_umi_file);
    }else{
        $UMI_RIGHT_number = 0;
    }
    #read-tag library
    my $LEFT_assemble_file = $every_dir.'/'.$LEFTSeqFile;
    my %sum_left_read;
    if (-s $LEFT_assemble_file){
        my $LEFT_assemble_handle = new IO::File "$LEFT_assemble_file" or die "$!";
        while (<$LEFT_assemble_handle>){
            chomp;
            my @line = split "\t",$_;
            if ($line[1] >= $reads_support_number){
                $sum_left_read{'good_num'} += 1;
                $sum_left_read{'good'} += $line[1];
            }else{
                $sum_left_read{'bad_num'} += 1;
                $sum_left_read{'bad'} += $line[1];
            }
        }
        $LEFT_assemble_handle -> close;
    }else{
        $sum_left_read{'good_num'} = 0;
        $sum_left_read{'bad_num'} = 0;
        $sum_left_read{'good'} = 0;
        $sum_left_read{'bad'} = 0;
    }

    my $RIGHT_assemble_file = $every_dir.'/'.$RIGHTSeqFile;
    my %sum_right_read;
    if (-s $RIGHT_assemble_file){
        my $RIGHT_assemble_handle = new IO::File "$RIGHT_assemble_file" or die "$!";
        while (<$RIGHT_assemble_handle>){
            chomp;
            my @line = split "\t",$_;
            if ($line[1] >= $reads_support_number){
                $sum_right_read{'good_num'} += 1;
                $sum_right_read{'good'} += $line[1];
            }else{
                $sum_right_read{'bad_num'} += 1;
                $sum_right_read{'bad'} += $line[1]
            }
        }
        $RIGHT_assemble_handle -> close;
    }else{
        $sum_right_read{'good_num'} = 0;
        $sum_right_read{'bad_num'} = 0;
        $sum_right_read{'good'} = 0;
        $sum_left_read{'bad'} = 0;
    }
    #Split16SFile
    my $split_file = $every_dir.'/'.$Split16SFile;
    my %split;
    if (-s $split_file){
        my $split_file_strings = gzip_support($split_file,'input',2);
        my $split_handle = new IO::File "$split_file_strings" or die "$!";
        while (my $head = <$split_handle>){
            my $seq = <$split_handle>;
            my $mid = <$split_handle>;
            my $qual = <$split_handle>;
            chomp $head;
            if ($head =~ /^@(.+)_\d+/){
                $split{$1} ++;
            }
        }
        $split_handle -> close;
    }else{
        warn "$every_dir\tnot $split_file or size eq 0!\n";
    }
    #Contig
    my $contig_file = $every_dir.'/'.$Contig;
    my %contig = ();
    my %contig2 = ();
    if (-s $contig_file){
        my $contig_handle = new IO::File "$contig_file" or die "$!";
        while (my $head = <$contig_handle>){
            my $seq = <$contig_handle>;
            chomp $seq;
            my $seq_length = length $seq;
            $contig{$seq_length} ++;
            $contig2{$seq_length} ++ if length $seq >= $long_contig_length;
        }
        $contig_handle -> close;
    }else{
        warn "$every_dir\tnot $contig_file or size eq 0!\n";
    }
    #cul N50
    my $contig_n50 = cul_n50(\%contig);
    #Clear Contig
    my $clear_contig_file = $every_dir.'/'.$ClearContig;
    my $clear_contig_number = 0;
    if (-s $clear_contig_file) {
        my $clear_contig_handle = new IO::File "$clear_contig_file" or die "$!";
        while (my $head = <$clear_contig_handle>){
            my $seq = <$clear_contig_handle>;
            chomp $seq;
            if (length $seq > 0){
                $clear_contig_number ++;
            }
        }
        $clear_contig_handle -> close;
    }
    #print
    my $split_number = scalar keys %split;
    my $contig_number = sum(values %contig);
    my $contig_gt_number = sum(values %contig2);
    print $output_handle "$dirname\t$UMI_LEFT_number\t$UMI_RIGHT_number\t$sum_left_read{'good_num'}\t$sum_left_read{'good'}\t$sum_left_read{'bad_num'}\t$sum_left_read{'bad'}\t$sum_right_read{'good_num'}\t$sum_right_read{'good'}\t$sum_right_read{'bad_num'}\t$sum_right_read{'bad'}\t$split_number\t$contig_number\t$contig_n50\t$contig_gt_number\t$clear_contig_number\n";
}
$output_handle -> close;

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
    my $input_file_handle = new IO::File "$input_file" or die "$!";
    while (<$input_file_handle>){
        $file_line_number ++;
    }
    $input_file_handle -> close;
    return $file_line_number;
}

sub cul_n50 {
    my $hash = shift;
    my %hash = %{$hash};
    my $total_number = sum(values %hash);
    my $n50_result = '';
    my $mid_number = '';
    my $mid_second = '';
    if ($total_number % 2 == 1){
        $mid_number = int($total_number / 2) + 1;
    }else{
        $mid_number = int($total_number / 2);
        $mid_second = int($total_number / 2) + 1;
    }
    my $now_index = 0;
    my @result;
    foreach my $every_length (sort {$a <=> $b} keys %hash){
        for (my $i = 1; $i <= $hash{$every_length}; $i ++){
            $now_index ++;
            if ($mid_second =~ /\d/){
                push @result,$every_length if $now_index == $mid_number or $now_index == $mid_second;
            }else{
                if ($now_index == $mid_number){
                    $n50_result = $every_length;
                }
            }
        }
    }
    if ($mid_second =~ /\d/){
        $n50_result = sum(@result) / 2;
    }
    return $n50_result;
}
