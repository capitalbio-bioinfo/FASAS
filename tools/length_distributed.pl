#!/usr/bin/perl
=head
    FilaName: Length_distributed.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Length_distributed.pl
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
    Input: Sequence file
    Output: Length statistcs file
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my ($input_file,$output_file,$type,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "type=s" => \$type,
    "help" => \$help,
) or die "Failed to Get Options!\n";

if ($help or $options_number == 0){
    print "Length Statistics Program of Sequence
Usage:
    perl $0 -i|--input_file [input file] -o|--output_file [output file] -t|--type [file format]
Options:
    input_file     strings    Input file [required]
    output_file    strings    Output file [required]
    type           strings    File format [fa|fq, required]
    help           switch     Print this message [off]

";
    exit;
}

die "LD: $input_file not read!\n" if ! -r $input_file;
die "LD: TYPE $type ERROR!\n" if $type !~ /fa/i and $type !~ /fq/i;

my $input_handle = new IO::File "$input_file" or die "$!";
my $output_handle = new IO::File ">$output_file" or die "$!";

my %dis;
if ($type =~ /fa/i or $type =~ /fasta/i){
    while (my $head = <$input_handle>){
        my $seq = <$input_handle>;
        chomp $seq;
        my $length = length $seq;
        $dis{$length} ++;
    }
}elsif ($type =~ /fq/i or $type =~ /fastq/i){
    while (my $head = <$input_handle>){
        my $seq = <$input_handle>;
        my $mid = <$input_handle>;
        my $qual = <$input_handle>;
        chomp $seq;
        my $length = length $seq;
        $dis{$length} ++;
    }
}
$input_handle -> close;

foreach my $every (sort { $a <=> $b } keys %dis){
    print $output_handle "$every\t$dis{$every}\n";
}
$output_handle -> close;
