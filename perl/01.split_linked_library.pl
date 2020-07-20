#!/usr/bin/perl
=head
    FileName: 01.Split_UMI_Library.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    01.Split_UMI_Library.pl : Part of the FASAS basic analysis script
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
    Input: Linked-Tag library all sequence
    Output: Correctly recognized sequences and bad sequences
=cut

=data_structure
FASTQ gzip or no gzip
R1:
    @.+
    [regular sequence][ATCGN]+
    \+
    FFFFF+
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;
use List::Util qw(min uniq);

my ($r1_file,$r2_file,$umi_library_r1_p,$max_distance,$log_file,$help);

my %dege = (
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

my %dege_out = (
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

#default options
$log_file = 'linked-tag_library.log';
$max_distance = 3;

my $options_number = scalar @ARGV;
GetOptions (
    "forward_file=s" => \$r1_file,
    "reverse_file=s" => \$r2_file,
    "umir1p=s" => \$umi_library_r1_p,
    "max_distance=i" => \$max_distance,
    "help" => \$help,
) or echo_error("SplitLinkedLibrary","Failed to get options");

if ($help or $options_number < 1){
    print "Linked-tag Library Split Program
Usage:
    perl $0 -f|--forward_file [file] -r|--reverse_file [file] -u|--umir1p [strings] -m|--max_distance [int]
Options:
    forward_file    strings    Linked-tag library forward sequence file         [none]
    reverse_file    strings    Linked-tag library reverse sequence file         [none]
    umir1p          strings    Linked-tag library front end regular sequence    [none]
    max_distance    int        The maximum distance of umir1p                   [3]
    help            switch     Show this help messages and exit                 [none]

";
    exit;
}

echo_error("SplitLinkedLibrary","$r1_file no exists!") if ! -e $r1_file;
echo_error("SplitLinkedLibrary","$r2_file no exists!") if ! -e $r2_file;
#echo_error("SplitLinkedLibrary","")

#thread ref sequence
chomp $umi_library_r1_p;
my $ur1p_length = length $umi_library_r1_p;
my @ur1p = split "",$umi_library_r1_p;
my @ur1p_out;
for (my $i = 0; $i <= $#ur1p; $i ++){
    if (exists $dege{$ur1p[$i]}){
        push @ur1p_out,$dege_out{$ur1p[$i]};
        $ur1p[$i] = $dege{$ur1p[$i]};
    }else{
        push @ur1p_out,$ur1p[$i];
    }
}
$umi_library_r1_p = join("",@ur1p);
my $umi_regular_r1_p = join("",@ur1p_out);

#get input handle
my $umi_r1_strings = gzip_support($r1_file,'input',2);
my $umi_r2_strings = gzip_support($r2_file,'input',2);
my $umi_r1_handle = new IO::File "$umi_r1_strings" or echo_error("SplitLinkedLibrary","$!");
my $umi_r2_handle = new IO::File "$umi_r2_strings" or echo_error("SplitLinkedLibrary","$!");

#get output handle
my $r1_out_strings = gzip_support('linked-tag_R1.fastq.gz','output',2);
my $r2_out_strings = gzip_support('linked-tag_R2.fastq.gz','output',2);
my $r1_out_handle = new IO::File "$r1_out_strings" or echo_error("SplitLinkedLibrary","$!");
my $r2_out_handle = new IO::File "$r2_out_strings" or echo_error("SplitLinkedLibrary","$!");

my $r1_fail_strings = gzip_support('Unmatch_linked-tag_R1.fastq.gz','output',2);
my $r2_fail_strings = gzip_support('Unmatch_linked-tag_R2.fastq.gz','output',2);
my $r1_fail_handle = new IO::File "$r1_fail_strings" or echo_error("SplitLinkedLibrary","$!");
my $r2_fail_handle = new IO::File "$r2_fail_strings" or echo_error("SplitLinkedLibrary","$!");
#print "$umi_r1_strings\t$umi_r2_strings\n";

my $log_handle = new IO::File ">$log_file" or echo_error("SplitLinkedLibrary","$!");

#defined stat variable
my $umi_pass_total = 0;
my $umi_total = 0;
my $pos_zero = 0;
my $pos_one = 0;
my $pos_two = 0;
my $correct_fail_num = 0;
my $correct_pass_num = 0;
my $correct_total_num = 0;

#go
while (my $head = <$umi_r1_handle>){
    my $seq = <$umi_r1_handle>;
    my $middle = <$umi_r1_handle>;
    my $qual = <$umi_r1_handle>;
    my $head2 = <$umi_r2_handle>;
    my $seq2 = <$umi_r2_handle>;
    my $middle2 = <$umi_r2_handle>;
    my $qual2 = <$umi_r2_handle>;
    chomp $seq;
    chomp $qual;
    my @SeqThread = SeqThread($seq);
    if ( $SeqThread[0] == 0 ){
        if (length $SeqThread[1] ne length $qual){
            my $qual_number = -(length $SeqThread[1]);
            $qual = substr $qual,$qual_number;
        }
        print $r1_out_handle "$head$SeqThread[1]","\n+\n","$qual\n";
        print $r2_out_handle "$head2$seq2$middle2$qual2";
        $umi_pass_total ++;
    }else{
        print $r1_fail_handle "$head$seq","\n+\n","$qual\n";
        print $r2_fail_handle "$head2$seq2$middle2$qual2";
    }
    $umi_total ++;
}
$umi_r1_handle -> close;
$umi_r2_handle -> close;
$r1_out_handle -> close;
$r2_out_handle -> close;
$r1_fail_handle -> close;
$r2_fail_handle -> close;

#print stat
my $successful_ratio = 100 * $umi_pass_total / $umi_total;
printf $log_handle "UMI Library Total Seq: %d, Split Successful Seq: %d,pr: %.2f%%.\n",$umi_total,$umi_pass_total,$successful_ratio;
my $ratio_zero = 100 * $pos_zero / $umi_pass_total;
my $ratio_one = 100 * $pos_one / $umi_pass_total;
my $ratio_two = 100 * $pos_two / $umi_pass_total;
printf $log_handle "UMI Library Pos Stat, Zero: %d,pr: %.2f%%. One: %d,pr: %.2f%%. Two: %d,pr: %.2f%%.\n",$pos_zero,$ratio_zero,$pos_one,$ratio_one,$pos_two,$ratio_two;
my $ratio_fail = 100 * $correct_fail_num / $correct_total_num;
my $ratio_pass = 100 * $correct_pass_num / $correct_total_num;
printf $log_handle "UMI Library Correct Stat, Correct total: %d,Correct Pass: %d,pr: %.2f%% ,Minimum distance gt $max_distance: %d,pr: %.2f%%.\n",$correct_total_num,$correct_pass_num,$ratio_pass,$correct_fail_num,$ratio_fail;

$log_handle -> close;

#sub
sub SeqThread {
    my $in_seq = shift;
    #limit postion: 1-3 bp
    #return 0 is success, return 1 is fail
    if ($in_seq =~ /^$umi_library_r1_p/){
        $pos_zero ++;
        return 0,$in_seq;
    }elsif ($in_seq =~ /\w($umi_library_r1_p.+)$/){
        my $ture_seq = $1;
        $pos_one ++;
        return 0,$ture_seq;
    }elsif ($in_seq =~ /\w\w($umi_library_r1_p.+)$/){
        $pos_two ++;
        my $ture_seq = $1;
        return 0,$ture_seq;
    }else{
        $correct_total_num ++;
        my %distance = ();
        for (my $i = 0;$i <= 2;$i ++){
            my $source_seq = substr $in_seq,$i,$ur1p_length;
            my $distance = HanMingDistance($source_seq);
            $distance{$i} = $distance if $distance <= $max_distance;
        }
        if (scalar keys %distance > 0){
            my $min_values = min(values %distance);
            my @min_number = grep { $distance{$_} == $min_values } keys %distance;
            if (scalar @min_number == 1){
                my $ture_seq_begin = $min_number[0] + $ur1p_length;
                my $ture_seq = substr $in_seq,$ture_seq_begin;
                $ture_seq = $umi_regular_r1_p.$ture_seq;
                $correct_pass_num ++;
                return 0,$ture_seq;
            }else{
                my $ture_seq_begin = 2 + $ur1p_length;
                my $ture_seq = substr $in_seq,$ture_seq_begin;
                $ture_seq = $umi_regular_r1_p.$ture_seq;
                $correct_pass_num ++;
                return 0,$ture_seq;
            }
        }else{
            $correct_fail_num ++;
            return 1,undef;
        }
    }
}

sub HanMingDistance {
    my $source_seq = shift;
    my $distance = 0;
    for (my $i = 0; $i <= $ur1p_length - 1; $i ++){
        $distance ++ if substr($source_seq,$i,1) !~ /${ur1p[$i]}/;
        return 100 if $distance > $max_distance;
    }
    return $distance;
}
