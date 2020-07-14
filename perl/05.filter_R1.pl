#!/usr/bin/perl
=head
    FileName: 05.Filter_R1.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    05.Filter_R1.pl : Part of the FASAS basic analysis script
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
    Input: Read-Tag library sequence
    Output: Result of filter
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

my ($input_r1,$input_r2,$output_r1,$output_r2,$type,$format,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_r1=s" => \$input_r1,
    "input_r2=s" => \$input_r2,
    "output_r1=s" => \$output_r1,
    "output_r2=s" => \$output_r2,
    "type=s" => \$type,
    "format=s" => \$format,
    "help" => \$help,
) or echo_error("FilterR1","Failed to get options");

if ($help or $options_number == 0){
    print "Assemble Library R1 Sequence Filter Program
Usage:
    perl $0 --input_r1 --input_r2 --output_r1 --output_r2 --format --append
Options:
    input_r1     strings    Input r1 file [required]
    input_r2     strings    Input r2 file [required]
    output_r1    strings    Output r1 file [required]
    output_r2    strings    Output r2 file [required]
    format       strings    Format of input file [fa|fq, required]
    help         switch     Print this message [off]

";
    exit;
}

my $type_word = 'FilterR1'.'-'.$type;

echo_error($type_word,"$input_r1 not be read!") if ! -r $input_r1;
echo_error($type_word,"$input_r2 not be read!") if ! -r $input_r2;

my $input_r1_strings = gzip_support($input_r1,'input',2);
my $input_r2_strings = gzip_support($input_r2,'input',2);
my $input_r1_handle = new IO::File "$input_r1_strings" or echo_error($type_word,"$!");
my $input_r2_handle = new IO::File "$input_r2_strings" or echo_error($type_word,"$!");

my $output_r1_strings = gzip_support($output_r1,'output',4);
my $output_r2_strings = gzip_support($output_r2,'output',4);
my $output_r1_handle = new IO::File "$output_r1_strings" or echo_error($type_word,"$!");
my $output_r2_handle = new IO::File "$output_r2_strings" or echo_error($type_word,"$!");

my $counter = 0;
my $total = 0;

if ($format =~ /fa/i or $format =~ /fasta/i){
    while (my $head = <$input_r1_handle>){
        my $seq = <$input_r1_handle>;
        my $head2 = <$input_r2_handle>;
        my $seq2 = <$input_r2_handle>;
        chomp $seq;
        $total ++;
        my $last_str = substr($seq,-20);
        my $G_num = () = $last_str =~ /G/ig;
        $counter ++ and next if $G_num >= 15;
        print $output_r1_handle "$head$seq\n";
        print $output_r2_handle "$head2$seq2";
    }
}elsif ($format =~ /fq/i or $format =~ /fastq/i){
    while (my $head = <$input_r1_handle>){
        my $seq = <$input_r1_handle>;
        my $mid = <$input_r1_handle>;
        my $qual = <$input_r1_handle>;
        my $head2 = <$input_r2_handle>;
        my $seq2 = <$input_r2_handle>;
        my $mid2 = <$input_r2_handle>;
        my $qual2 = <$input_r2_handle>;
        chomp $seq;
        $total ++;
        my $last_str = substr($seq,-20);
        my $G_num = () = $last_str =~ /G/ig;
        $counter ++ and next if $G_num >= 15;
        print $output_r1_handle "$head$seq\n$mid$qual";
        print $output_r2_handle "$head2$seq2$mid2$qual2";
    }
}else{
    echo_error("GCRatio","format is $format, Unknown");
}

$input_r1_handle -> close;
$input_r2_handle -> close;
$output_r1_handle -> close;
$output_r2_handle -> close;

echo_error("$type_word","Total sequence is zero") if $total == 0;
my $pr = sprintf "%0.4f%%",$counter / $total * 100;

echo_info("$type_word","Delete sequence $counter, pr: $pr");
