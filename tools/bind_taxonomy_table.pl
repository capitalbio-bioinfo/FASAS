#!/usr/bin/perl
=head
    FileName: Bind_TaxonomyTable.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Bind_TaxonomyTable.pl
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
    Input: FASAS project folder
    Output: Annotation file at a certain rank
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;
use File::Basename;

my ($input_dir,$type,$output,$number_output,$filter_word,$help);

my $options_number = scalar @ARGV;
GetOptions(
    'input_dir=s' => \$input_dir,
    'type=s' => \$type,
    'output=s' => \$output,
    'number' => \$number_output,
    'filter_word=s' => \$filter_word,
    'help' => \$help,
);

if ($help or $options_number == 0){
    print "Bind the Same Level Annotation Files in a 16S-FASAS Project
Usages:
    perl $0 -i|--input_dir [project folder] -t|--type [anno level] -o|--output [output file]
Options:
    input_dir    path       Input path, a 16S-FASAS Project [required]
    type         strings    Annotation Level, ex: Genus or Species [required]
    output       file       Output file [required]
    number       switch     Output Number or Decimal [default: Decimal]
    filter_word  strings    Filter samples using regular expressions[default: no]
    help         switch     Print this message [off]
";
    exit;
}

my $now_output_file = $input_dir.'/'.$output;
unlink $now_output_file if -e $now_output_file;
my $common_name = $type.'.tsv';
my $file = `find $input_dir |grep -i $common_name`;
my @file = split "\r?\n",$file;

#filter
if (defined $filter_word){
    @file = grep { $_ !~ /$filter_word/ } @file;
}

my %tax;
my %total;
my %species;
foreach my $every_file (@file){
    my $strings = '_'.$common_name;
    my $sample_name = basename($every_file,$strings);
    #$tax{$sample_name}{}
    my $sample_handle = new IO::File "$every_file" or die "$!";
    while (<$sample_handle>){
        chomp;
        next if $_ =~ /^Tax/;
        my @line = split "\t",$_;
        $tax{$sample_name}{$line[0]} = $line[1];
        $total{$sample_name} += $line[1];
        $species{$line[0]} = 1 if ! exists $species{$line[0]};
    }
    $sample_handle -> close;
}

#thread tax
if (! $number_output){
    foreach my $every_sample (keys %total){
        foreach my $every_species (keys %{$tax{$every_sample}}){
            $tax{$every_sample}{$every_species} = $tax{$every_sample}{$every_species} / $total{$every_sample};
        }
    }

    my $output_handle = new IO::File ">$output" or die "$!";
    print $output_handle "Tax";
    foreach my $every_sample (sort { $a cmp $b } keys %tax){
        print $output_handle "\t$every_sample";
    }
    print $output_handle "\n";
    foreach my $every_species (keys %species){
        print $output_handle "$every_species";
        foreach my $every_sample (sort { $a cmp $b } keys %tax){
            if (exists $tax{$every_sample}{$every_species}){
                printf $output_handle "\t%.8f",$tax{$every_sample}{$every_species};
            }else{
                print $output_handle "\t0";
            }
        }
        print $output_handle "\n";
    }
    $output_handle -> close;
}else{
    my $output_handle = new IO::File ">$output" or die "$!";
    print $output_handle "Tax";
    foreach my $every_sample (sort { $a cmp $b } keys %tax){
        print $output_handle "\t$every_sample";
    }
    print $output_handle "\n";
    foreach my $every_species (keys %species){
        print $output_handle "$every_species";
        foreach my $every_sample (sort { $a cmp $b } keys %tax){
            if (exists $tax{$every_sample}{$every_species}){
                print $output_handle "\t$tax{$every_sample}{$every_species}";
            }else{
                print $output_handle "\t0";
            }
        }
        print $output_handle "\n";
    }
    $output_handle -> close;
}
