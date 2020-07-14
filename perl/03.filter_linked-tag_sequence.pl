#!/usr/bin/perl
=head
    FileName: 03.Filter_UMISeq.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    03.Filter_UMISeq.pl : Part of the FASAS basic analysis script
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
    Input: Assemble Linked-Tag library sequence
    Output: Filter Linked-Tag library sequence
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

echo_info("FilterUMISeq","BEGIN");

my ($upper_limit,$input_file,$output_file,$help);

my $option_number = scalar @ARGV;
GetOptions (
    "upper_limit=i" => \$upper_limit,
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "help" => \$help,
) or die "Failed to get options!\n";

if ($help or $option_number == 0){
    print "UMI Library Sequence Length Filter Program
Usage:
    perl $0 --upper_limit [int] --input_file [input_file] --output_file [output_file]
Options:
    upper_limit   int     The threshold of the sequence length [required]
    input_file    path    input file, FASTQ format [required]
    output_file   path    output file, FASTA format [required]
    help          switch  get help information [off]

";
    exit;
}

echo_error("FilterUMISeq","$input_file cannot be read!") if ! -r $input_file;

my $input_handle = new IO::File "$input_file" or echo_error("FilterUMISeq","$!");
my $output_handle = new IO::File ">$output_file" or echo_error("FilterUMISeq","$!");

my $counter = 0;
my $total = 0;
while (my $head = <$input_handle>){
    $head =~ s/^@/>/;
    my $seq = <$input_handle>;
    my $mid = <$input_handle>;
    my $qual = <$input_handle>;
    chomp $seq;
    $total ++;
    $counter ++ and next if length $seq > $upper_limit;
    print $output_handle "$head$seq\n";
}
echo_error("FilterUMISeq","Sequence number is zero") if $total == 0;
my $pr = sprintf "%0.4f%%",$counter / $total * 100;

$input_handle -> close;
$output_handle -> close;

echo_info("FilterUMISeq","Filter Sequence $counter, pr: $pr");
echo_info("FilterUMISeq","END");
