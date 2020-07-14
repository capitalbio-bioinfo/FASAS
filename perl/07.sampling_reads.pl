#!/usr/bin/perl
=head
    FileName: 07.PE_Read_Sampling.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    07.PE_Read_Sampling.pl : Part of the FASAS basic analysis script
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
    Input: All 16s rRNA sequence
    Output: 16s rRNA sampling data
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

my ($input_r1,$input_r2,$output_r1,$output_r2,$sampling_number,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_r1=s" => \$input_r1,
    "input_r2=s" => \$input_r2,
    "output_r1=s" => \$output_r1,
    "output_r2=s" => \$output_r2,
    "sampling_number=i" => \$sampling_number,
    "help" => \$help,
) or echo_error("Failed to get options");

if ($help or $options_number == 0){
    print "Pair End Sequencing Data Sampling Program
Usages:
    perl $0 --input_r1 [r1] --input_r2 [r2] --output_r1 [r1] --output_r2 [r2] --sampling_number [int]
Options:
    input_r1           strings    Input r1 file [required]
    input_r2           strings    Input r2 file [required]
    output_r1          strings    Output r1 file [required]
    output_r2          strings    Output r2 file [required]
    sampling_number    int        Sampling number [required]
    help               switch     Print this massage [off]

";
    exit;
}

echo_error("PEReadSampling","$input_r1 cannot be read!") if ! -r $input_r1;
echo_error("PEReadSampling","$input_r2 cannot be read!") if ! -r $input_r2;

#get file col number and create rand hash
my $total_sequence_number = 0;
my $number_strings = gzip_support($input_r1,'input',2);
my $number_handle = new IO::File "$number_strings" or echo_error("PEReadSampling","$!");

while (<$number_handle>){
    $total_sequence_number ++;
}
$number_handle -> close;

$total_sequence_number = int($total_sequence_number / 4);
my %rand;

if ($sampling_number >= $total_sequence_number){
    for (my $i = 1; $i <= $total_sequence_number; $i ++){
        $rand{$i} = 1;
    }
}else{
    srand();
    for (my $i = 1; $i <= $sampling_number; $i ++){
        RERAND:
        my $rand_num = int(rand($total_sequence_number));
        goto RERAND if exists $rand{$rand_num};
        $rand{$rand_num} = 1;
    }
}
#sampling
my $input_r1_strings = gzip_support($input_r1,'input',2);
my $input_r2_strings = gzip_support($input_r2,'input',2);
my $input_r1_handle = new IO::File "$input_r1_strings" or echo_error("PEReadSampling","$!");
my $input_r2_handle = new IO::File "$input_r2_strings" or echo_error("PEReadSampling","$!");

my $output_r1_strings = gzip_support($output_r1,'output',4);
my $output_r2_strings = gzip_support($output_r2,'output',4);
my $output_r1_handle = new IO::File "$output_r1_strings" or echo_error("PEReadSampling","$!");
my $output_r2_handle = new IO::File "$output_r2_strings" or echo_error("PEReadSampling","$!");

#read fasta
my $number = 1;
my $print_number = 0;
while (my $r1_head = <$input_r1_handle>) {
    my $r1_seq = <$input_r1_handle>;
    my $r2_head = <$input_r2_handle>;
    my $r2_seq = <$input_r2_handle>;
    next if $r1_head !~ /^[@>]/ or $r2_head !~ /^[@>]/;
    $r1_head =~ s/^@/>/;
    $r2_head =~ s/^@/>/;
    if (exists $rand{$number}){
        print $output_r1_handle "$r1_head$r1_seq";
        print $output_r2_handle "$r2_head$r2_seq";
        $print_number ++;
    }
    $number ++;
}
$input_r1_handle -> close;
$input_r2_handle -> close;
$output_r1_handle -> close;
$output_r2_handle -> close;

echo_warn("PEReadSampling","The number of sequences is less than $sampling_number, the result of sequence coverage check may be inaccurate\n") if $print_number < $sampling_number * 0.5;
