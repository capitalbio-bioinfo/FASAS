#!/usr/bin/perl
=head
    FileName: 06.StatAssembleRead.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    06.StatAssembleRead.pl : Part of the FASAS basic analysis script
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
    Input: Filtered Read-Tag library sequence file
    Output: Read-Tag library statistics data
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

my ($reverse_fastq,$type,$output_number,$output_seqname,$help);

my $opt_num = scalar @ARGV;
GetOptions (
    "reverse_fastq=s" => \$reverse_fastq,
    "type=s" => \$type,
    "output_number=s" => \$output_number,
    "output_seqname=s" => \$output_seqname,
    "help" => \$help,
) or echo_error("StatAssemble","Failed to get options");

if ($help or $opt_num == 0){
    print "Stat Assemble Library Read Program
Usage:
    perl $0 -r [input_r2_fastq] -t [input_type]
Options:
    r    strings    input file, FASTQ format [required]
    t    strings    type of input file [required]
";
    exit;
}


echo_error('StatAssemble',"Unknown type: $type") if $type !~ /l/i and $type !~ /r/i;
echo_error('StatAssemble',"Input Data is emtpy") if ! defined $reverse_fastq;
echo_error('StatAssemble',"$reverse_fastq not found") if ! -e $reverse_fastq;

my $type_word;
my $output_tsv;
if ($type =~ /l/i){
    $type_word = 'StatAssemble-LEFT';
}else{
    $type_word = 'StatAssemble-RIGHT';
}

#input
my $input_strings = gzip_support($reverse_fastq,'input',2);
my $input_handle = new IO::File "$input_strings" or echo_error("$type_word",$!);

#output
##seqname
my $seqname_handle = new IO::File ">$output_seqname" or echo_error("$type_word",$!);
##number
my $number_handle = new IO::File ">$output_number" or echo_error("$type_word",$!);

my %hash;
my %name;
while (<$input_handle>){
    chomp;
    my @header = split " ",$_;
    my $seq = <$input_handle>; chomp $seq;
    my $middle = <$input_handle>; chomp $middle;
    my $qual = <$input_handle>; chomp $middle;
    if ($hash{$header[1]}){
        $hash{$header[1]} ++;
        $name{$header[1]} .= "\t".$header[0];
    }else{
        $hash{$header[1]} = 1;
        $name{$header[1]} = $header[0];
    }
}
$input_handle -> close;

foreach my $every_key (sort {$hash{$b} <=> $hash{$a}} keys %hash){
    print $number_handle $every_key,"\t",$hash{$every_key},"\n";
}

foreach my $every_key (keys %name){
    print $seqname_handle $every_key,"\t",$name{$every_key},"\n";
}
$seqname_handle -> close;
$number_handle -> close;
