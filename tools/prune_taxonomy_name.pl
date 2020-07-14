#!/usr/bin/perl
=head
    FileName: Cut_LevelAnnotate.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Cut_LevelAnnotate.pl
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
    Input: Annotation table
    Output: Changed annotation table
=cut

=example
    ex. one line of input data:
        Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides stercoris	114
    ex. one line of output data:
        Bacteroides stercoris	114
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my ($input_tax_file,$taxonomy_level,$output_tax_file,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_tax_file=s" => \$input_tax_file,
    "output_tax_file=s" => \$output_tax_file,
    "taxonomy_level=s" => \$taxonomy_level,
    "help" => \$help,
) or die "Failed to Get Options!\n";

if ($help or $options_number == 0){
    print "Annotation Information Simplification Program
Usage:
    perl $0 -i [input taxonomy file] -o [output taxonomy file]
Options:
    input_tax_file     file       input taxonomy file [required]
    output_tax_file    file       output taxonomy file [required]
    taxonomy_level     strings    taxonomy level [required]
    help               switch     Print this message [off]

";
    exit;
}

my $input_handle = new IO::File "$input_tax_file" or die "$!";
my $output_handle = new IO::File ">$output_tax_file" or die "$!";

my $header = <$input_handle>;
print $output_handle "$header";

while (<$input_handle>){
    chomp;
    my @line = split "\t",$_;
    $line[0] = newTax($line[0]);
    print $output_handle join("\t",@line),"\n";
}

$input_handle -> close;
$output_handle -> close;

sub newTax {
    my $tax = shift;
    my @tax = split ";",$tax;
    my $last_tax = $tax[-1];
    my $before_tax = '';
    for (my $i = $#tax; $i >= 0; $i --){
        if ($tax[$i] =~ /__/ or $tax[$i] =~ /PAC\d+/ or $tax[$i] =~ /no taxonomy/){
            next;
        }else{
            $before_tax = $tax[$i];
            last;
        }
    }
    print "$last_tax\t$before_tax\n";
    if ($before_tax eq $last_tax){
        return $before_tax;
    }else{
        my $new_tax = $last_tax.'_'.$before_tax;
        return $new_tax;
    }
}
