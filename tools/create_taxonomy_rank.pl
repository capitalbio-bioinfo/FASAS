#!/usr/bin/perl
=head
    FileName: Create_TaxonomyRank.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Create_TaxonomyRank.pl
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
    Output: Rank statistcs table for all sample
=cut

use strict;
use warnings;
use apackages;
use IO::File;
use Getopt::Long;
use File::Basename qw(basename);
use File::Find;
use Cwd qw(abs_path);
use List::Util qw(uniq);

my ($search_folder,$output_file,$rank_array,$help);

#default
our @rank_array = qw(kingdom phylum class order family genus species);
$rank_array = '2,3,4,5,6,7';

my $option_number = scalar @ARGV;
GetOptions (
    "search_folder=s" => \$search_folder,
    "output_file=s" => \$output_file,
    "rank_array=s" => \$rank_array,
    "help" => \$help,
) or echo_error("CreateTaxonomyRank","Failed to get options!");

#check options
if ($help or $option_number < 1){
    print "
Usage:
    perl $0 --search_folder|-s [16S-FASAS project folder] --output_file|-o [output file]
Options:
    search_folder    strings    
    output_file      strings    
    rank_array       strings    [2,3,4,5,6,7]
    help             switch     print this messages [off]
";
    exit;
}

#thread rank array
if ($rank_array){
    echo_error("CreateTaxonomyRank","Other character: $1") if $rank_array =~ /([^\d,])/;
    my @filter_array = split ",",$rank_array;
    my $other_strings = grep { $_ !~ /^\d$/ } @filter_array;
    if ($other_strings){
        echo_error("CreateTaxonomyRank","Option value ERROR of rank_array, the correct ex: -r 1,2,3");
    }
    my %sure_array;
    foreach (@filter_array){
        $sure_array{$rank_array[$_ - 1]} = 1;
    }
    @rank_array = grep { defined $sure_array{$_} } @rank_array;
}

our @all_file = '';
find(\&find_condition, $search_folder);

#filter file
@all_file = grep { filterfile($_) } @all_file;
@all_file = map { abs_path($_) } @all_file;
@all_file = uniq(@all_file);

#stat rank
my %result;
foreach my $every_file (@all_file){
    my $file_name = basename($every_file);
    if ($file_name =~ /(.+)_([a-zA-Z]+)\.tsv/){
        my $smaple_name = $1;
        my $file_rank = $2;
        my $rank_number = wcfile($every_file);
        $result{$smaple_name}{$file_rank} = $rank_number;
    }else{
        echo_error("CreateTaxonomyRank","$every_file ???");
    }
}

#output
my $output_handle = new IO::File ">$output_file" or echo_error("CreateTaxonomyRank","Failed to open output handle!");
print $output_handle "Sample\t";
print $output_handle join("\t",@rank_array),"\n";

foreach my $every_sample (sort {$a cmp $b} keys %result){
    print $output_handle "$every_sample";
    foreach my $every_rank (@rank_array){
        echo_error("CreateTaxonomyRank","Sample: $every_sample, Rank: $every_rank, annotation file not found!") if $result{$every_sample}{$every_rank} !~ /\d+/;
        print $output_handle "\t$result{$every_sample}{$every_rank}";
    }
    print $output_handle "\n";
}
$output_handle -> close;

sub find_condition {
    push @all_file,$File::Find::name;
}

sub wcfile {
    my $address = shift;
    my $input_handle = new IO::File "$address" or echo_error("CreateTaxonomyRank","Failed to open handle!");
    my $row_number = -1;
    while (<$input_handle>){
        next if $_ =~ /^#/;
        $row_number ++;
    }
    $input_handle -> close;
    return $row_number;
}

sub filterfile {
    my $file = shift;
    foreach my $every_rank (@rank_array){
        if ( $file =~ /_$every_rank\.tsv/ ){
            return 1;
        }
    }
    return 0;
}
