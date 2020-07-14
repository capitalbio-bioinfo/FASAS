#!/usr/bin/perl
=head
    FileName: 00.Analysis_Library.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    00.Analysis_Library.pl : Part of the FASAS basic analysis script
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
    Input: Library info file
    Output: Strings of library info
=cut

use strict;
use warnings;
use Cwd qw(abs_path);
use Switch;
use IO::File;
use Getopt::Long;
use File::Basename;
use apackages;

my ($input_file, $help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "help" => \$help,
) or echo_error("AnalysisLibrary","Failed to get options");

if ($help or $options_number < 1){
    print "Library Structure Analysis Program
Usage:
    perl $0 -i [input_file] -u [umi_length]
Options:
    input_file   strings     A text file describing library design    [none]
    help         switch      Show this help message and exit          [off]

";
    exit;
}

echo_error("AnalysisLibrary","Library Structure File $input_file not found!") if ! -e $input_file;

my %library_info = (
    'LinkedLibraryR1P' => '',
    'LinkedLibraryR2P' => '',
    'LinkedLibraryR1M' => '',
    'LinkedLibraryR2M' => '',
    'LinkedLibraryR1A' => '',
    'LinkedLibraryR2A' => '',
    'ReadLibraryLeftP' => '',
    'ReadLibraryRightP' => '',
    'ReadLibraryLeftA' => '',
    'ReadLibraryRightA' => '',
    'UMILength' => '',
);

#find umi length
my $umi_length;
my $input_handle = new IO::File "$input_file" or echo_error("AnalysisLibrary","$!");
while (<$input_handle>){
    chomp;
    next if $_ =~ /^#/;
    my @line = split " ",$_;
    my @umi = $line[1] =~ /N+/g;
    foreach my $every_umi (@umi){
        next if length $every_umi < 7;
        if (! $umi_length){
            $umi_length = length $every_umi;
        }else{
            echo_error("AnalysisLibrary","AnalysisLibraryInfo: The length of the UMI is not unique!
Please check library info file: $input_file") if $umi_length != length $every_umi;
        }
    }
}
$input_handle -> close;

echo_error("AnalysisLibrary","AnalysisLibraryInfo: Not Found UMI!\nPlease check library info file: $input_file") if $umi_length eq '';
$library_info{'UMILength'} = $umi_length;

#read info txt
$input_handle = new IO::File "$input_file" or echo_error("AnalysisLibrary","$!");
while (<$input_handle>){
    chomp;
    next if $_ =~ /^#/;
    my @line = split " ",$_;
    switch ($line[0]) {
        case { $line[0] =~ /R1/i } { ThreadR1Library($line[1]) }
        case { $line[0] =~ /R2/i } { ThreadR2Library($line[1]) }
        case { $line[0] =~ /Left/i } { ThreadReadLibrary('left',$line[1]) }
        case { $line[0] =~ /Right/i } { ThreadReadLibrary('right',$line[1]) }
        else { die "ERROR: unrecognized library type! R1 or R2\n$_" }
    };
}
$input_handle -> close;

#print array
my @array_queue = qw(UMILength LinkedLibraryR1P LinkedLibraryR1M LinkedLibraryR1A LinkedLibraryR2P LinkedLibraryR2M LinkedLibraryR2A ReadLibraryLeftP ReadLibraryLeftA ReadLibraryRightP ReadLibraryRightA);

my @array = ();
foreach my $every_array (@array_queue){
    push @array,$library_info{$every_array};
}

print join(" ",@array);

sub ThreadR1Library {
    my $input_seq = shift;
    #dege N not split, need length
    my $umi_n = 'N' x $umi_length;
    my @regular_seq = split "$umi_n",$input_seq;
    #return $regular_seq[0],$regular_seq[1],$regular_seq[2],
    if (scalar @regular_seq == 3){
        #UMI library
        $library_info{'LinkedLibraryR1P'} = $regular_seq[0];
        $library_info{'LinkedLibraryR1M'} = $regular_seq[1];
        $library_info{'LinkedLibraryR1A'} = $regular_seq[2];
    }elsif (scalar @regular_seq == 2){
        $library_info{'ReadLibraryR1P'} = $regular_seq[0];
        $library_info{'ReadLibraryR1A'} = $regular_seq[1];
    }else{
        echo_error("AnalysisLibrary","ERROR: Regular Seq Number is Abnormal!");
    }
}

sub ThreadR2Library {
    my $input_seq= shift;
    #dege N not split, need length
    my $umi_n = 'N' x $umi_length;
    my @regular_seq = split "$umi_n",$input_seq;
    if (scalar @regular_seq == 3){
        #Linked-tag Library
        $library_info{'LinkedLibraryR2P'} = $regular_seq[0];
        $library_info{'LinkedLibraryR2M'} = $regular_seq[1];
        $library_info{'LinkedLibraryR2A'} = $regular_seq[2];
    }elsif (scalar @regular_seq == 2){
        #Read-tag Library
        $library_info{'ReadLibraryR2P'} = $regular_seq[0];
        $library_info{'ReadLibraryR2A'} = $regular_seq[1];
    }else{
        echo_error("AnalysisLibrary","ERROR: Regular Seq Number is Abnormal!");
    }
}

sub ThreadReadLibrary {
    my $seq_type = shift;
    my $input_seq = shift;
    my $umi_n = 'N' x $umi_length;
    my @regular_seq = split "$umi_n",$input_seq;
    if (scalar @regular_seq == 2){
        if ($seq_type eq 'left'){
            $library_info{'ReadLibraryLeftP'} = $regular_seq[0];
            $library_info{'ReadLibraryLeftA'} = $regular_seq[1];
        }elsif ($seq_type eq 'right'){
            $library_info{'ReadLibraryRightP'} = $regular_seq[0];
            $library_info{'ReadLibraryRightA'} = $regular_seq[1];
        }else{
            echo_error("AnalysisLibrary","ERROR: Seq type is $seq_type, expect left or right!");
        }
    }else{
        echo_error("AnalysisLibrary","ERROR: Regular seq number is abnormal!");
    }
}
