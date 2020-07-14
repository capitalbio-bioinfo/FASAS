#!/usr/bin/perl
=head
    FileName: 08.Stat_ContigLength.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    08.Stat_ContigLength.pl : Part of the FASAS basic analysis script
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
    Input: Nogap contigs
    Output: Statistcs data of contigs length
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

my ($min_value,$max_value,$input_file,$output_file,$help);

my $option_number = scalar @ARGV;
GetOptions (
    "min_value=i" => \$min_value,
    "max_value=i" => \$max_value,
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "help" => \$help,
) or echo_error("StatContigLength","Failed to get options");

if ($help or $option_number == 0){
    print "Stat Contig Length Program
Usage:
    perl $0 --min_value [int] --max_value [int] --input_file [input_file] --output_file [output_file]
Options:
    min_value     int     minimum value of the statistics file
    max_value     int     maximum value of the statistics file
    input_file    path    input file, FASTQ format [required]
    output_file   path    output file, FASTA format [required]
    help          switch  get help information [off]

";
    exit;
}

echo_error("StatContigLength","$input_file cannot be read!") if ! -r $input_file;

my $input_handle = new IO::File "$input_file" or echo_error("StatContigLength","$!");
my $output_handle = new IO::File ">$output_file" or echo_error("StatContigLength","$!");

my %hash;
while (my $head = <$input_handle>){
    my $seq = <$input_handle>;
    chomp $seq;
    my $length = length $seq;
    $hash{$length} ++ if $length >= $min_value and $length <= $max_value;
}
$input_handle -> close;

foreach my $every_pos (sort { $a <=> $b } keys %hash){
    print $output_handle "$every_pos\t$hash{$every_pos}\n";
}
$output_handle -> close;
