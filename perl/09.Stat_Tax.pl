#!/usr/bin/perl
=head
    FileName: 09.Stat_Tax.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    09.Stat_Tax.pl : Part of the FASAS basic analysis script
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
    Input: Primer annotation file
    Output: Annotation table for each rank
=cut

use strict;
use warnings;
use IO::File;
use Switch;
use apackages;
use Getopt::Long;

my ($sample_name,$input_file,$output_dir,$precent,$help);

my $option_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "output_dir=s" => \$output_dir,
    "sample_name=s" => \$sample_name,
    "precent" => \$precent,
    "help" => \$help,
) or die "Failed to Get Options!";

if ($help or $option_number == 0){
    print "Get Annotate Infomations Program
Usage:
    perl $0 -i [input_file] -o [output_dir] -s [sample_name]
Options:
    input_file    path      input file, result file form 09.Annotate_Tax.pl script [required]
    output_dir    path      output dir [required]
    sample_name   strings   sample name [required]
    precent       switch    conversion number is decimal [off]
    help          switch    get help information [off]

";
    exit;
}

#check options
die "GetAnnoInfo: Input File $input_file cannot be read!\n" if ! -r $input_file;
die "GetAnnoInfo: Sample Name \"$sample_name\" cannot be read!\n" if $sample_name eq '';

my %tax;
my @level = qw(kingdom phylum class order family genus species);
my $input_handle = new IO::File "$input_file" or die "$!";
mkdir $output_dir if ! -d $output_dir;

my $sum_number = 0;
#input format: contig_name	id	tax
while (<$input_handle>){
    chomp;
    my @line = split "\t",$_;
    if (scalar @line < 3){
        echo_warn("StatTax","$_ looks like only two columns, omit!");
        next;
    }
    if ($line[2] eq 'unclassified'){
        my $unclassified_cycle = 7 - 1 - 1;
        for (my $i = 0; $i <= $unclassified_cycle; $i ++){
            $line[2] .= ';unclassified';
        }
    }
    my @tax = split ";",$line[2];
    if (scalar @tax < 7 and $line[2] ne 'unclassified'){
        echo_warn("StatTax","$_ classification less than 7 levels, omit!");
        next;
    }elsif (scalar @tax > 7){
        echo_warn("StatTax","$_ classification more than 7 levels, omit!");
        next;
    }
    $sum_number ++;
    for (my $i=0; $i<=$#tax; $i++){
        my $str = join(";",@tax[0..$i]);
        switch ($i) {
            case 0 { $tax{'kingdom'}{$str} ++ }
            case 1 { $tax{'phylum'}{$str} ++ }
            case 2 { $tax{'class'}{$str} ++ }
            case 3 { $tax{'order'}{$str} ++ }
            case 4 { $tax{'family'}{$str} ++ }
            case 5 { $tax{'genus'}{$str} ++ }
            case 6 { $tax{'species'}{$str} ++ }
        }
    }
}
$input_handle -> close;

##filter only one contigs
#foreach my $every_rank (keys %tax){
#    foreach my $every_tax (keys %{$tax{$every_rank}}){
#        if($tax{$every_rank}{$every_tax} == 1){
#            delete $tax{$every_rank}{$every_tax};
#        }
#    }
#}

#prepare handle
my %output;
foreach my $every_level (@level){
    my $output_file = $output_dir.'/'.$sample_name.'_'.$every_level.'.tsv';
    $output{$every_level} = new IO::File ">$output_file" or die "$!";
    print { $output{$every_level} } "Tax\t$sample_name\n";
}

#primer print
if ($precent){
    foreach my $every_key (keys %tax){
        foreach my $second_key (keys %{$tax{$every_key}}){
            printf { $output{$every_key} } "%s\t%.8f\n",$second_key,$tax{$every_key}{$second_key} / $sum_number;
        }
    }
}else{
    foreach my $every_key (keys %tax){
        foreach my $second_key (keys %{$tax{$every_key}}){
            print { $output{$every_key} } "$second_key\t$tax{$every_key}{$second_key}\n";
        }
    }
}

foreach my $every_handle (keys %output){
    $output{$every_handle} -> close;
}
