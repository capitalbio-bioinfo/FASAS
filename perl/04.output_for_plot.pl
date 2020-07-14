#!/usr/bin/perl
=head
    FileName: 04.OutputForPlot.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    04.OutputForPlot.pl : Part of the FASAS basic analysis script
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
    Input: UMI table
    Output: UMI statistical data
    Why: Rscript takes up too much memory when dealing with very large tables
=cut

=algorithm
Read UMI Table
Keep the first two columns
Filter data by the second column gt 2
=cut

use strict;
use warnings;
use apackages;
use IO::File;
use Getopt::Long;

my ($input_file,$output_file,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "h" => \$help,
) or echo_error("OutputForPlot","Failed to get options");

if ($help or $options_number < 2){
    print "Output For Rscript Plot Program
Usage:
    perl $0 -i StatUMI_LtoR.tsv -o [output file]
Options:
    input_file     strings    input file, umi table file in step4
    output_file    strings    output file, left umi or right umi file for Rscript
    help           switch     output this message
";
    exit;
}

echo_error("OutputForPlot","$input_file no exists!") if ! -e $input_file;
echo_error("OutputForPlot","output_file is NA!") if ! defined $output_file;

my $input_handle = new IO::File "$input_file" or echo_error("OutputForPlot","$!");
my $output_handle =  new IO::File ">$output_file" or echo_error("OutputForPlot","$!");

while (<$input_handle>){
    chomp;
    my @line = split "\t",$_;
    #if lt 2, filter line
    next if $line[1] <= 2;
    print $output_handle "$line[0]\t$line[1]\n";
}

$input_handle -> close;
$output_handle -> close;

__END__
