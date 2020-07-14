#!/usr/bin/perl
=head
    FileName: 00.Analysis_parameter.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    00.Analysis_Parameter.pl : Part of the FASAS basic analysis script
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
    Input: ConfigFile
    Output: Opyions Array
=cut

=configFileStrings
    SampleName
    WorkFolder
    LibraryInfo
    LinkedLibraryR1
    LinkedLibraryR2
    ReadLibraryR1
    ReadLibraryR2
    CoverageDatabase
    ReferenceFasta
    ReferenceTaxonomy
    AssembleProgram
    IlluminaAdapter
    ContigLength
=cut

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use apackages;

my ($input_file,$help);
my %parameter;

my $option_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "help" => \$help,
) or echo_error("AnalysisParameter","Failed to get options");

if ($help or $option_number == 0){
    print "FASAS Parameter Analysis Program
Usage:
    perl $0 --input_file [config file]
Options:
    input_file    strings    A FASAS ConfigFile                 [none]
    help          switch     Show this help message and exit    [off]
";
    exit;
}

#Read ConfigFile
my $input_handle = new IO::File "$input_file" or echo_error("AnalysisParameter","$!");
my $line_number = 0;
while (<$input_handle>){
    $line_number ++;
    next if $_ =~ /^#/;
    next if $_ =~ /^\s+?$/;
    chomp;
    if ($_ =~ /(.+)\s+=\s+(.+)/ or $_ =~ /(.+)\s?=\s?(.+)/){
        $parameter{$1} = $2;
    }else{
        echo_error("AnalysisParameter","Line $line_number regular expression processing failed\n$_");
    }
}
$input_handle -> close;

#Process parameters in order
my @array = qw(SampleName WorkFolder LibraryInfo LinkedLibraryR1 LinkedLibraryR2 ReadLibraryR1 ReadLibraryR2 CoverageDatabase ReferenceFasta ReferenceTaxonomy AssembleProgram IlluminaAdapter ContigLength);
my @parameter;
foreach my $every_value (@array){
    if (exists $parameter{$every_value}){
        push @parameter,$parameter{$every_value};
        delete $parameter{$every_value};
    }else{
        echo_error("AnalysisParameter","Parameter $every_value is missing");
    }
}

#Process unclear parameters
if (scalar keys %parameter > 0){
    foreach my $every_key (keys %parameter){
        echo_warn("AnalysisParameter","Parameter $every_key don't understand!");
    }
}

#Output to STDOUT
print join(' ',@parameter);
