#!/usr/bin/perl
=head
    FileName: 03.Correct_UMISeq.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    03.Correct_UMISeq.pl : Part of the FASAS basic analysis script
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
    Input: ConfigFile
    Output: Opyions Array
=cut

=note
UMI Library Sequence Correction Program base on Simple Dynamic Alignment
For 16S-FASAS Design
Input: UMI Library Fasta File
InputFileStatus:
    --input_fasta Fasta File should be quality controlled, included QC, Merge and Length Filter
Output: The Correct Fasta File
OutputFileStatus:
    --output_file Primer Output File
LOG: Common File, TEXT Format
LOGFileStatus:
    -low_similarity_log Record the Reasons and Alignmnet that do not meet the Conditional Sequence
    -umi_toolong_log Record the Sequence that UMI Length gt 14
=cut

=testlog
20190430: version 1.0.1, 566.172K read, 40 threads, 5 min 17 second, ratio of primer result: 99.81%, UMI too long ratio: 0.1574%
=cut

use strict;
use warnings;
use threads;
use Fcntl qw(:flock);
use Switch;
use apackages;
use POSIX qw(strftime);
use IO::File;
use File::Path qw(make_path remove_tree);
use Getopt::Long;
use List::Util qw(max);

my %dege = (
    'A' => 'A',
    'T' => 'T',
    'C' => 'C',
    'G' => 'G',
    'R' => '[AG]',
    'Y' => '[CT]',
    'M' => '[AC]',
    'K' => '[GT]',
    'S' => '[GC]',
    'W' => '[AT]',
    'H' => '[ATC]',
    'B' => '[GTC]',
    'V' => '[GAC]',
    'D' => '[GAT]',
    'N' => '[ATCG]',
);

#In order not to output degenerate bases
my %trandege = (
    'R' => 'A',
    'Y' => 'C',
    'M' => 'A',
    'K' => 'G',
    'S' => 'G',
    'W' => 'A',
    'H' => 'A',
    'B' => 'G',
    'V' => 'G',
    'D' => 'G',
    'N' => 'A',
);

#defined default options
my ($input_fasta,$output_file,$umi_toolong_log,$low_similarity_log,$help,$ref_seq,$ref_length,$ur1p,$ur1m,$ur1a);
my $umi_length = 14;
my $max_threads = 40;
my $match_score = 2;
my $N_score = 1;
my $mismatch_score = -2;
my $gap_score = -8;

#get options
my $option_number = scalar @ARGV;
GetOptions (
    "input_fasta=s" => \$input_fasta,
    "output_file=s" => \$output_file,
    "umi_toolong_log=s" => \$umi_toolong_log,
    "low_similarity_log=s" => \$low_similarity_log,
    "ur1p=s" => \$ur1p,
    "ur1m=s" => \$ur1m,
    "ur1a=s" => \$ur1a,
    "umi_length=i" => \$umi_length,
    "threads=i" => \$max_threads,
    "match=i" => \$match_score,
    "nmatch=i" => \$N_score,
    "mismatch=i" => \$mismatch_score,
    "gap=i" => \$gap_score,
    "help" => \$help,
) or echo_error("CorrectUMISeq","Failed to get options");

if ($help or $option_number == 0){
    print "UMI Library Sequence Correction Program
    input_fasta        path       input file, FASTA format [required]
    output_file        path       output file, FASTA format [required]
    umi_toolong_log    path       sequence filtered because UMI is too long, FASTA format [required]
    low_similarity_log           path       log file, TEXT format [required]
    ur1p               sequence   regular sequence in UMI Library [required]
    ur1m               sequence   regular sequence in UMI Library [required]
    ur1a               sequence   regular sequence in UMI Library [required]
    umi_length         int        umi length in UMI library [14]
    threds             int        max threads number [40]
    match              int        base match score [2]
    nmatch             int        UMI match score [1]
    mismatch           int        base mismatch penalty [-2]
    gap                int        gap penalty [-8]
    help               switch     get help information [off]
";
    exit 0;
}

#check options
echo_error("CorrectUMISeq","Parameter input_fasta is required") if $input_fasta eq '';
echo_error("CorrectUMISeq","Parameter output_file is required") if $output_file eq '';
echo_error("CorrectUMISeq","Parameter umi_toolong_log is required") if $umi_toolong_log eq '';
echo_error("CorrectUMISeq","Parameter low_similarity_log is required") if $low_similarity_log eq '';
echo_error("CorrectUMISeq","Input: $input_fasta not found!") if ! -r $input_fasta;

echo_info("CorrectUMISeq","Create reference sequence and reference regular expression");
#create ref_seq
chomp ($ur1p,$ur1m,$ur1a);
my $umi_n = 'N' x $umi_length;
$ref_seq = $ur1p.$umi_n.$ur1m.$umi_n.$ur1a;

#create regular expression
my $umi_lower_limit = $umi_length - 1;
my $umi_upper_limit = $umi_length + 1;
my $ref_regular = $ur1p.'(.{'.$umi_lower_limit.','.$umi_upper_limit.'})'.$ur1m.'(.{'.$umi_lower_limit.','.$umi_upper_limit.'})'.$ur1a;
my @ref_regular = split "",$ref_regular;
$ref_regular = join('',map { exists $dege{$_} ? $dege{$_} : $_ } @ref_regular);
undef @ref_regular;

#get ref info
$ref_length = length $ref_seq;
my $ur1p_length = length $ur1p;
my $ur1m_length = length $ur1m;
my $ur1a_length = length $ur1a;
my $ur1p_out = TranSeq($ur1p);
my $ur1m_out = TranSeq($ur1m);
my $ur1a_out = TranSeq($ur1a);

#get umi position
my $umi_one_begin = $ur1p_length;
my $umi_two_begin = $ur1p_length + $umi_length + $ur1m_length;

#create temp dir
my $local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
my $tmpdir = '/tmp/Correct_UMISeq_'.$local_time.'_TmpDir';
make_path($tmpdir) ? echo_info("CorrectUMISeq","Create TMP Dir $tmpdir") : echo_error("CorrectUMISeq","Failed to Create TMP Dir");

my $available_cpu_threads = CheckCPU();
$max_threads = $available_cpu_threads if $available_cpu_threads < $max_threads;
echo_info("CorrectUMISeq","Correct Begin, Threads Number is $max_threads");
#threads
my @threads;
my %data;
my $thread_number = 0;
my $read_index = 1;
#open input handle
my $input_fasta_hadnle = new IO::File "$input_fasta" or echo_error("CorrectUMISeq","$!");
while (my $header = <$input_fasta_hadnle>){
    my $seq = <$input_fasta_hadnle>;
    chomp $header;
    chomp $seq;
    $data{$header} = $seq;
    if ($read_index == 5000){
        $read_index = 1;
        $thread_number ++;
        if ($thread_number > $max_threads){
            my $multiple = int($thread_number/$max_threads);
            $thread_number = $thread_number - $max_threads * $multiple;
        }
        push @threads,threads -> create(\&ThreadsDNA,\%data,$thread_number,$tmpdir);
        %data = ();
        if (scalar @threads == $max_threads){
            my $every_thread = shift @threads;
            $every_thread -> join();
        }
    }
    $read_index ++;
}
$input_fasta_hadnle -> close;

#analysis last data
echo_info("CorrectUMISeq","The program is analyzing the final data");
if (scalar keys %data > 0){
    $thread_number ++;
    if ($thread_number > $max_threads){
        my $multiple = int($thread_number/$max_threads);
        $thread_number = $thread_number - $max_threads * $multiple;
    }
    push @threads,threads -> create(\&ThreadsDNA,\%data,$thread_number,$tmpdir);
}

#join remaining threads
while (my $every_thread = shift @threads){
    $every_thread -> join();
}

echo_info("CorrectUMISeq","Merge result file and log file into working directory");
#cat result and log file
my $umi_tmp_file = $tmpdir.'/thread_*.fasta';
my $umi_tmp_file2 = $tmpdir.'/thread_2_*.log';
my $log_tmp_file = $tmpdir.'/thread_*.log';
my $cat_status = system("cat $umi_tmp_file > $output_file");
$cat_status == 0 ? echo_info("CorrectUMISeq","Cat Result to $output_file") : echo_error("CorrectUMISeq","Failed to cat TMP Result!");
$cat_status = system("cat $umi_tmp_file2 > $umi_toolong_log");
$cat_status == 0 ? echo_info("CorrectUMISeq","Cat Length LOG to $umi_toolong_log") : echo_error("CorrectUMISeq","Failed to cat TMP Length LOG!");
$cat_status = system("cat $log_tmp_file > $low_similarity_log");
$cat_status == 0 ? echo_info("CorrectUMISeq","Cat LOG to $umi_toolong_log") : echo_error("CorrectUMISeq","Failed to cat TMP LOG!");

#remove temp dir
echo_info("CorrectUMISeq","Remover tmp folder");
remove_tree($tmpdir) ? echo_info("CorrectUMISeq","Remove TMP Dir $tmpdir") :  echo_error("CorrectUMISeq","Remove TMP Dir ERROR");

echo_info("CorrectUMISeq","END");

sub ThreadsDNA {
    my %data = %{ shift @_ };
    my $thread_num = shift;
    my $tmpdir = shift;
    my $output_tmp_file = $tmpdir.'/thread_'.$thread_num.'.fasta';
    my $output_tmp_file2 = $tmpdir.'/thread_2_'.$thread_num.'.log';
    my $log_tmp_file = $tmpdir.'/thread_'.$thread_num.'.log';
    my $output_tmp_handle = new IO::File ">>$output_tmp_file" or echo_error("$output_tmp_file:$!");
    flock $output_tmp_handle,LOCK_EX;
    my $length_log_handle = new IO::File ">>$output_tmp_file2" or echo_error("$output_tmp_file2:$!");
    flock $length_log_handle,LOCK_EX;
    my $log_tmp_handle = new IO::File ">>$log_tmp_file" or echo_error("$log_tmp_file:$!");
    flock $log_tmp_handle,LOCK_EX;
    foreach my $every_header (keys %data){
        #step1 match possible sequence
        if ($data{$every_header} =~ /$ref_regular/){
                printf $output_tmp_handle "%s\n%s\n",$every_header,$data{$every_header};
                next;
        }
        #step2 correction
        my ($match_num,$algin_result_one,$algin_result_two) = AlignSeq($data{$every_header});
        if ($match_num > $ref_length * 0.8){
            my @analysis_result = AnalysisDA($algin_result_one,$algin_result_two);
            my $umi_one_seq = substr($data{$every_header},$analysis_result[0],$analysis_result[1]);
            my $umi_two_seq = substr($data{$every_header},$analysis_result[2],$analysis_result[3]);
            if (length $umi_one_seq <= $umi_length and length $umi_two_seq <= $umi_length){
                #chioce $output_tmp_handle
                printf $output_tmp_handle "%s\n%s%s%s%s%s\n",$every_header,$ur1p_out,$umi_one_seq,$ur1m_out,$umi_two_seq,$ur1a_out;
            }else{
                #gt 14, delete seq
                #chioce $length_log_handle
                printf $length_log_handle "%s\n%s\n",$every_header,$data{$every_header};
            }
        }else{
            print $log_tmp_handle "$every_header: Match Number is $match_num < 80% of the Ref Sequence Length\n";
            print $log_tmp_handle "Alignmnet Result:\n";
            print $log_tmp_handle "$algin_result_one\n";
            print $log_tmp_handle "$algin_result_two\n";
        }
    }
    flock $output_tmp_handle,LOCK_UN;
    flock $length_log_handle,LOCK_UN;
    flock $log_tmp_handle,LOCK_UN;
    $output_tmp_handle -> close;
    $length_log_handle -> close;
    $log_tmp_handle -> close;
}

sub AnalysisDA {
    my $align_ref = shift;
    my $align_query = shift;
    my $umi_one_now_begin = $umi_one_begin;
    my $umi_two_now_begin = $umi_two_begin;
    my $umi_one_length = $umi_length;
    my $umi_two_length = $umi_length;
    my $ref_num = 1;
    my $query_num = 1;
    if (length $align_ref == length $align_query){
        for (my $i = 0; $i < $ref_length; $i ++){
            if (substr($align_ref,$i,1) eq '-'){
                switch ($ref_num) {
                    case {$ref_num < $umi_one_begin} {$umi_one_now_begin ++ && $umi_two_now_begin ++}
                    case {$ref_num == $umi_one_begin} {$umi_two_now_begin ++ && $umi_one_length ++}
                    case {$ref_num < $umi_two_begin} {$umi_two_now_begin ++}
                    case {$ref_num == $umi_two_begin} {$umi_one_length ++}
                }
            }else{
                $ref_num ++;
            }
        }
        for (my $j = 0; $j < $ref_length; $j ++){
            if (substr($align_query,$j,1) eq '-'){
                switch ($query_num) {
                    case {$query_num <= $umi_one_begin} {$umi_one_now_begin -- && $umi_two_now_begin --}
                    case {$query_num <= $umi_two_begin} {$umi_two_now_begin --}
                }
            }else{
                $query_num ++;
            }
        }
    }else{
        echo_error("WHAT FOUK IS THIS!!!\n$align_ref\n$align_query");
    }
    return ($umi_one_now_begin,$umi_one_length,$umi_two_now_begin,$umi_two_length);
}

sub AlignSeq {
    my $input_seq = shift;
    my %score;
    my $input_length = length $input_seq;
    for (my $i=0; $i <= $ref_length; $i ++){
        for (my $j=0; $j <= $input_length; $j++){
            $score{$j}{$i} = 0;
        }
    }
    #initialization
    ##this step is very fast
    for (my $i = 0;$i <= $ref_length; $i ++){
        $score{0}{$i} = $i * $gap_score;
    }
    for (my $j = 0;$j <= $input_length; $j ++){
        $score{$j}{0} = $j * $gap_score;
    }
    #calculated score matrix
    ##speed limit step one
    for (my $i = 1;$i <= $ref_length; $i ++){
        for (my $j = 1;$j <= $input_length; $j ++){
            $score{$j}{$i} = max($score{$j-1}{$i-1} + AlignBase(substr($ref_seq,$i-1,1),substr($input_seq,$j-1,1)),$score{$j-1}{$i} + $gap_score,$score{$j}{$i-1} + $gap_score);
        }
    }
    #traceback
    ##speed limit step two
    my $align_ref = '';
    my $align_query = '';
    my ($i,$j) = (length $ref_seq,length $input_seq);
    my $match_num = 0;
    while ($i > 1 or $j > 1){
        if ($i == 0){
            $align_ref = '-'.$align_ref;
            $j -= 1;
        }elsif ($j == 0){
            $align_query = '-'.$align_query;
            $i -= 1;
        }elsif ($score{$j}{$i} == $score{$j-1}{$i} + $gap_score){
            $align_ref = '-'.$align_ref;
            $align_query = substr($input_seq,$j-1,1).$align_query;
            $j -= 1;
        }elsif ($score{$j}{$i} == $score{$j}{$i-1} + $gap_score) {
            $align_query = '-'.$align_query;
            $align_ref = substr($ref_seq,$i-1,1).$align_ref;
            $i -= 1;
        }else {
            $match_num ++;
            $align_ref = substr($ref_seq,$i-1,1).$align_ref;
            $align_query = substr($input_seq,$j-1,1).$align_query;
            $j -= 1;
            $i -= 1;
        }
    }
    return ($match_num,$align_ref,$align_query);
}

sub AlignBase {
    my $ref_base = shift;
    my $two_base = shift;
    return $N_score if $ref_base eq 'N';
    $two_base =~ /$dege{$ref_base}/ ? return $match_score : return $mismatch_score;
}

sub TranSeq {
    my $input_seq = shift;
    my @input_seq = split "",$input_seq;
    for (my $i = 0; $i <= $#input_seq; $i ++){
        $input_seq[$i] = $trandege{$input_seq[$i]} if exists $trandege{$input_seq[$i]};
    }
    $input_seq = join("",@input_seq);
    return $input_seq;
}
