#!/usr/bin/perl
=head
    FileName: Base_balance.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Base_balance.pl : Part of the FASAS basic analysis script
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
    Output: Sequence base statistcs data
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;
use apackages;

my ($input_file,$output_file,$format,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "format=s" => \$format,
    "help" => \$help,
) or die "Failed to Get Options!\n";

if ($help or $options_number == 0){
    print "Base Statistics Program of Sequence
Usage:
    perl $0 -i|--input_file [input file] -o|--output_file [output file] -f|--format [input file format]
Options:
    input_file      path       Input file [required]
    output_file     path       Output file [default: STDOUT]
    format          strings    Input file format [fa|fq, required]
    help            switch     Print this message [off]

";
    exit;
}

my $output_handle;
if ($output_file){
    $output_handle = new IO::File ">$output_file" or die "$!";
}else{
    $output_handle = <STDOUT>;
}

my %hash;
my %total;
my @array = qw(A T C G);
my $fasta_handle = new IO::File "$input_file" or die "$!";
while (my($id, $decs, $seq) = read_complex_fast($fasta_handle)){
    $total{$id} = length $seq;
    foreach my $every_array (@array){
        $hash{$id}{$every_array} = () = $seq =~ /$every_array/ig;
    }
}
$fasta_handle -> close;

foreach my $every (keys %hash){
    print $output_handle $every;
    foreach my $every_sub (sort { $a cmp $b } keys %{$hash{$every}}){
        next if $every_sub !~ /[ATCG]/i;
        printf $output_handle "\t%s\t%.2f",$every_sub,$hash{$every}{$every_sub} / $total{$every};
    }
    print $output_handle "\n";
}
$output_handle -> close;
