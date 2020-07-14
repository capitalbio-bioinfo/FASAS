#!/usr/bin/perl
=head
    FileName: Filter_Megablast_result.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Filter_Megablast_result.pl
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
    Input: Taxonomy table
    Output: Filtered taxonomy table
=cut

use strict;
use warnings;
use IO::File;
use List::Util qw(sum);
use Getopt::Long;

my ($input_file,$output_file,$filter_value,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "filter_value=s" => \$filter_value,
    "help" => \$help,
) or die "$!";

my $input_handle = new IO::File "$input_file" or die "$!";
my $output_handle = new IO::File ">$output_file" or die "$!";

#read taxonomy table, N sample
my @file = <$input_handle>;
$input_handle -> close;

#get sample line and filter data
@file = grep { $_ !~ /^#/ } @file;
my $header = shift @file;
chomp $header;

my @result;
foreach my $every_tax (@file){
    chomp $every_tax;
    my @line = split "\t",$every_tax;
    my $tax_name = shift @line;
    my $tax_abundance = sum(@line);
    if ($tax_abundance > $filter_value){
        my $result_line = join("\t",$tax_name,@line);
        push @result,$result_line;
    }else{
        next;
    }
}

print $output_handle "$header\n";
foreach my $every_result (@result){
    print $output_handle "$every_result\n";
}
$output_handle -> close;
