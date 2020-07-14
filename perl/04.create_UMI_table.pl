#!/usr/bin/perl
=head
    FileName: 04.Create_UMITable.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    04.Create_UMITable.pl : Part of the FASAS basic analysis script
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
    Input: Clear Linked-Tag library sequence
    Output: A UMI statistical table
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;

my ($input_fasta,$left_table,$right_table,$gzip,$log_file,$help);
my ($ur1p,$ur1m,$ur1a,$umi_length);

my %dege = (
    'R' => '[AG]',
    'Y' => '[CT]',
    'M' => '[AC]',
    'K' => '[GT]',
    'S' => '[GC]',
    'W' => '[AT]',
    'H' => '[ATC]',
    'B' => '[GTC]',
    'V' => '[GAC]',
    'D' => '[GAT]',
    'N' => '[ATCG]',
);

my $option_number = scalar @ARGV;
GetOptions (
    "f=s" => \$input_fasta,
    "l=s" => \$left_table,
    "r=s" => \$right_table,
    "ur1p=s" => \$ur1p,
    "ur1m=s" => \$ur1m,
    "ur1a=s" => \$ur1a,
    "umi_length=i" => \$umi_length,
    "log_file=s" => \$log_file,
    "gzip" => \$gzip,
    "help" => \$help,
) or echo_error("CreateUMITable","Failed to get options");

if ($help or $option_number == 0){
    print "Create UMI Table Program
Usage:
    perl $0 -f [input_fasta] -l [left_table] -r [right_table] --ur1p [ur1p] --ur1m [ur1m] --ur1a [ur1a] --umi_length [number] --log_file [log_file]
Options:
    f           strings     input file, FASTA format [required]
    l           strings     output file, TSV format [required]
    r           strings     output file, TSV format [required]
    ur1p        sequence    regular sequence in UMI Library [required]
    ur1m        sequence    regular sequence in UMI Library [required]
    ur1a        sequence    regular sequence in UMI Library [required]
    umi_length  int         umi length in UMI library [14]
    log_file    strings     record the sequence regular match failed [required]
    gzip        switch      input_fasta is gzip format? [off]
    h           switch      get help information [off]
";
    exit;
}

#check options
echo_error("CreateUMITable","$input_fasta not found!") if ! -e $input_fasta;
echo_error("CreateUMITable","UR1P is \"$ur1p\", too short!") if length $ur1p < 4;
echo_error("CreateUMITable","UR1M is \"$ur1m\", too short!") if length $ur1m < 4;
echo_error("CreateUMITable","UR1A is \"$ur1a\", too short!") if length $ur1a < 4;
echo_error("CreateUMITable","UR1P is \"$ur1p\", has some illegal characters!") if $ur1p =~ /\s/;
echo_error("CreateUMITable","UR1M is \"$ur1m\", has some illegal characters!") if $ur1m =~ /\s/;
echo_error("CreateUMITable","UR1A is \"$ur1a\", has some illegal characters!") if $ur1a =~ /\s/;
echo_error("CreateUMITable","umi length is $umi_length, this is abnormal!") if $umi_length < 4 or $umi_length > 20;

my $log_handle = new IO::File ">$log_file" or echo_error("CreateUMITable","$!");
print $log_handle "#This file records the sequence in which the regular match failed.
#Most regular expression matching failures are due to UMI length error\n";
my $left_table_handle = new IO::File ">$left_table" or die "$!";
my $right_table_handle = new IO::File ">$right_table" or die "$!";

#thread ref seq
$ur1p = ThreadRef($ur1p);
$ur1m = ThreadRef($ur1m);
$ur1a = ThreadRef($ur1a);

#left_umi	right_umi	left <-> right
my @read_seq = ReadSequence($input_fasta);

my %umi_left_link = %{$read_seq[0]};
my %umi_right_link = %{$read_seq[1]};

#print table
#left
my %left_total_num;
foreach my $every_left (keys %umi_left_link){
    foreach my $every_right (keys %{$umi_left_link{$every_left}}){
        $left_total_num{$every_left} += $umi_left_link{$every_left}{$every_right};
    }
}
foreach my $every_left (sort { $left_total_num{$b} <=> $left_total_num{$a} } keys %umi_left_link){
    print $left_table_handle "$every_left\t$left_total_num{$every_left}";
    foreach my $every_right (sort { $umi_left_link{$every_left}{$b} <=> $umi_left_link{$every_left}{$a} } keys %{$umi_left_link{$every_left}}){
        print $left_table_handle "\t$every_right\t$umi_left_link{$every_left}{$every_right}";
    }
    print $left_table_handle "\n";
}
#right
my %right_total_num;
foreach my $every_right (keys %umi_right_link){
    foreach my $every_left (keys %{$umi_right_link{$every_right}}){
        $right_total_num{$every_right} += $umi_right_link{$every_right}{$every_left};
    }
}
foreach my $every_right (sort { $right_total_num{$b} <=> $right_total_num{$a} } keys %umi_right_link){
    print $right_table_handle "$every_right\t$right_total_num{$every_right}";
    foreach my $every_left (sort { $umi_right_link{$every_right}{$b} <=> $umi_right_link{$every_right}{$a} } keys %{$umi_right_link{$every_right}}){
        print $right_table_handle "\t$every_left\t$umi_right_link{$every_right}{$every_left}";
    }
    print $right_table_handle "\n";
}

$left_table_handle -> close;
$right_table_handle -> close;
$log_handle -> close;

sub ReadSequence {
    my $input_file = shift;
    my %umi_link_left;
    my %umi_link_right;
    my $input_handle;
    if ($gzip){
        $input_handle = new IO::File "gzip -dc $input_file|" or echo_error("CreateUMITable","$!");
    }else{
        $input_handle = new IO::File "$input_file" or echo_error("CreateUMITable","$!");
    }
    while (my $seq_header = <$input_handle>){
        my $seq_sequence = <$input_handle>;
        chomp $seq_sequence;
        if ($seq_sequence =~ /$ur1p([ATCGN]{$umi_length})$ur1m([ATCGN]{$umi_length})$ur1a/){
            $umi_link_left{$1}{$2} ++;
            $umi_link_right{$2}{$1} ++;
        }else{
            print $log_handle "$seq_sequence\n";
        }
    }
    $input_handle -> close;
    return \%umi_link_left,\%umi_link_right;
}

sub ThreadRef {
    my $input_seq = shift;
    chomp $input_seq;
    my @ref_array = split "",$input_seq;
    for (my $i = 0;$i <= $#ref_array;$i ++){
        $ref_array[$i] = $dege{$ref_array[$i]} if exists $dege{$ref_array[$i]};
    }
    $input_seq = join("",@ref_array);
    return $input_seq;
}
