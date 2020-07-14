#!/usr/bin/perl
=head
    FileName: GC_ratio.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    GC_ratio.pl
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
    Output: GC ratio of sequence statistcs data
=cut

use strict;
use warnings;
use threads;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use IO::File;
use Getopt::Long;
use File::Basename;
use Fcntl qw(:flock);
use File::Path qw(make_path remove_tree);


my ($input_file,$output_file,$threads_number,$format,$append,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_file=s" => \$input_file,
    "output_file=s" => \$output_file,
    "threads_number=i" => \$threads_number,
    "format=s" => \$format,
    "append" => \$append,
    "help" => \$help,
) or die "Failed to Get Options!\n";

if ($help or $options_number == 0){
    print "GC Ratio of Sequence Multithreading Statistics Program
Usage:
    perl $0 -i [input file] -o [output file] -t [threads number] -f [file format] -a [append]
Options:
    input_file      strings    Input file [required]
    output_file     strings    Output file [required]
    threads_number  int        Threads number [required]
    format          strings    Input file format [fa|fq, required]
    append          switch     Appends the result to the file [off]
    help            switch     Print this message [off]

";
    exit;
}

die "GCRatio: $input_file not be read!\n" if ! -r $input_file;
die "GCRatio: threads_number is $threads_number, too large" if $threads_number > 40;
my $file_row = `wc -l $input_file|awk '{print \$1}'`;
chomp $file_row;
my $d_row = $file_row / $threads_number;
die "GCRatio: InputFile Row.Number is $file_row, Threads.Number is $threads_number, Not divisible!\n" if int($d_row) != $d_row;

#append
if (-e $output_file and $append){
    print "GCRatio: The result appedded to $output_file.\n";
}elsif (-e $output_file and ! $append){
    print "GCRatio: Delete existing result files.\n";
    unlink $output_file;
}

my $local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
my $tmpdir = '/tmp/GCRatio_'.$local_time.'_TmpDir';
make_path($tmpdir) ? print "GCRatio: Create TMP Dir $tmpdir\n" : die "ERROR: Failed to Create TMP Dir \n";
GetSplitFile($tmpdir);
my $input_file_cwd = abs_path($input_file);
my $split_stat = system("bash $tmpdir/split.sh $tmpdir $input_file_cwd $threads_number");
die "GCRatio: split.sh return $split_stat!\n" if ! $split_stat == 0;

my $basename = basename($input_file);
my $out_basename = basename($output_file);
my $tmp_output_file = $tmpdir.'/'.$out_basename;
my @threads;
for (my $i = 0; $i < $threads_number; $i ++){
    my $file;
    if ($i < 10){
        $file = $tmpdir.'/'.$basename.'0'.$i;
    }else{
        $file = $tmpdir.'/'.$basename.$i;
    }
    die "GCRatio: $file not found!i is $i\n" if ! -e $file;
    push @threads,threads -> create(\&ThreadSeq,$file,$tmp_output_file,$format,$i);
}

while (my $every_thread = shift @threads){
    $every_thread -> join();
}

if ($append){
    system("cat $tmp_output_file >> $output_file");
}else{
    system("mv $tmp_output_file $output_file");
}
remove_tree($tmpdir) ? print "Remove TMP Dir $tmpdir\n" :  die "GCRatio: Remove TMP Dir ERROR\n";

sub ThreadSeq {
    my $seq_file = shift;
    my $output_result = shift;
    my $seq_format = shift;
    my %hash;
    my %total;
    my $input_handle = new IO::File "$seq_file" or die "$!";
    if ($seq_format =~ /fa/i or $seq_format =~ /fasta/i){
        while (my $head = <$input_handle>){
            my $seq = <$input_handle>;
            chomp $head;
            chomp $seq;
            $head =~ s/^>//;
            $total{$head} = length $seq;
            $hash{$head} = () = $seq =~ /[GC]/ig;
        }
    }elsif ($seq_format =~ /fq/i or $seq_format =~ /fastq/i){
        while (my $head = <$input_handle>){
            my $seq = <$input_handle>;
            my $mid = <$input_handle>;
            my $qual = <$input_handle>;
            chomp $head;
            chomp $seq;
            $head =~ s/^@//;
            $total{$head} = length $seq;
            $hash{$head} = () = $seq =~ /[GC]/ig;
        }
    }else{
        die "GCRatio: format is $format, Unknown"
    }
    my $output_handle = new IO::File ">>$output_result" or die "$!";
    flock $output_handle,LOCK_EX;
    foreach my $every (keys %hash){
        printf $output_handle "%s\t%.4f\n",$every,$hash{$every} / $total{$every};
    }
    flock $output_handle,LOCK_UN;
    $output_handle -> close;
}

sub GetSplitFile {
    my $tmpdir = shift;
    my $sh_file = $tmpdir.'/split.sh';
    my $sh_handle = new IO::File ">$sh_file" or die "$!";
    print $sh_handle '#!/bin/bash
cd $1
BaseName=$(basename $2)
ln -s $2 ./${BaseName}
RowNumber=$(wc -l ${BaseName} | cut -d" " -f 1)
EveryRowNumber=$((${RowNumber}/$3))
split -d -a 2 -l ${EveryRowNumber} ${BaseName} ${BaseName}
';
    $sh_handle -> close;
}
