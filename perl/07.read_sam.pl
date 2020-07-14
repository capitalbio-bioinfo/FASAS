#!/usr/bin/perl
=head
    FileName: 07.Read_Sam.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    07.Read_Sam.pl : Part of the FASAS basic analysis script
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
    Input: Bowtie2 SAM file
    Output: 16s rRNA Coverage
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

my ($input_path,$output_file,$type,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_path=s" => \$input_path,
    "output_file=s" => \$output_file,
    "type=s" => \$type,
    "help" => \$help,
) or echo_error("ReadSAM","Failed to Get Options!");

if ($help or $options_number < 1){
    print "Read SAM Script
Usage:
    perl $0 -i|--input_path [strings] -o|--output_file [strings] -t|--type [strings]
Options:
    input_path     strings    SAM file folder    [none]
    output_file    strings    
    type           strings    
    help           switch     
";
    exit;
}

die echo_error("ReadSAM","$input_path cannot be read!") if ! -e $input_path;
die echo_error("ReadSAM","Unknown sequence type!") if $type !~ /r1/i and $type !~ /r2/i;

my @input_file = @{ FindInputFile($input_path,$type) };

#initialization hash
#In order to prevent the value of NA from appearing in R
my %pos_stat = ();
for (my $i = 0; $i <= 1600; $i ++){
    $pos_stat{$i} = 0;
}

foreach my $every_sam (@input_file){
    ThreadSam($every_sam);
}

my $output_handle = new IO::File ">$output_file" or echo_error("ReadSAM","$!");
#create cov file
foreach my $every_pos (keys %pos_stat){
    print $output_handle "$every_pos , $pos_stat{$every_pos}\n";
}
$output_handle -> close;

sub ThreadSam {
    my $input_sam = shift;
    my $input_handle = new IO::File "$input_sam" or echo_error("ReadSAM","$!");
    while (<$input_handle>) {
        chomp;
        next if $_ =~ /^@/;
        my @sam_line = split "\t",$_;
        next if $sam_line[1] == 4;
        my $pos = $sam_line[3];
        my @CIGAR = $sam_line[5] =~ /\d+[A-Za-z]/g;
        my ($length,$status) = $CIGAR[0] =~ /(\d+)([A-Za-z])/;
        $pos = $pos - $length if $status =~ /S/i;
        for (my $i = 0; $i < length $sam_line[9]; $i ++){
            my $now_pos = $pos + $i;
            $pos_stat{$now_pos} ++;
        }
    }
    $input_handle -> close;
}

sub FindInputFile {
    my $path = shift;
    my $t = shift;
    my @file = ();
    opendir "DIR","$path" or echo_error("ReadSAM","Failed to open folder");
    while (my $every_file = readdir(DIR)){
        push @file,$every_file if $every_file =~ /sam$/i and $every_file =~ /$t/i;
    }
    closedir DIR;
    echo_error("ReadSAM","Did not find the required file!") if scalar @file < 1;
    return \@file;
}
