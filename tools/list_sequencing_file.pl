#!/usr/bin/perl
=head
    FileName: List_SequencingFile.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    List_SequencingFile.pl
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
    Output: Sample info file
=cut

=input_format
    #No Header in the Input File
    #The input file has only three columns
    #The first three columns are 'Sequencing ID', 'Sample ID' and 'Library Type'
    #Sequencing ID: String common to R1 and R2 of the sequencing file in front of the file name
    #Sample ID: sample name, one read-tag library and one linked-tag library share a sample name
    #Library Type: Read-tag library or Linked-tag library
=cut

use strict;
use warnings;
use Cwd qw(abs_path);
use IO::File;
use Getopt::Long;
use Data::Dumper;

my ($file_address,$help,$input_file,$split_word,$output_file);
$split_word = '';

my $options_number = scalar @ARGV;
GetOptions (
    "file_address=s" => \$file_address,
    "input_file=s" => \$input_file,
    "split_word=s" => \$split_word,
    "output_file=s" => \$output_file,
    "help" => \$help,
) or die "$!";

if ($help or $options_number < 1){
    print "Sequence File Address Acquisition Program
Usage:
    perl $0 -i [sample_info.tsv] -f [path] -s [split word]
Options:
    file_address    path        Directory for storing sequence files
    input_file      path        Sample info table, format: SequencingID,SampleID,LibraryType
    output_file     path        Output file
    split_word      strings     Symbol for splitting the sample name
    help            switch      Output this message
";
    exit;
}

$file_address = abs_path($file_address);
die "Address Folder no exists!" if ! -d $file_address;
die "$input_file not found!" if ! -e $input_file;

my $input_handle = new IO::File "$input_file" or die "$!";
my $output_handle = new IO::File ">$output_file" or die "$!";

#read sample information and obtain sequencing file
my %hash; # sequencing id => [ sample, librarytype, R1, R2 ]
while (<$input_handle>){
    chomp;
    next if $_ =~ /^#/;
    my @line = split "\t",$_;
    my $name = `find $file_address |grep $line[0] | sort`;
    my @name = split "\r?\n",$name;
    if (defined $split_word){
        @name = grep { $_ =~ /$line[0]$split_word/ } @name;
    }else{
        my $a_w = '[^0-9A-Za-z]';
        @name = grep { $_ =~ /$line[0]$a_w/ } @name;
    }
    die "$_\nFailed to found sequencing file!\n" if scalar @name ne 2;
    unshift @name,$line[2];
    unshift @name,$line[1];
    $hash{$line[0]} = [@name];
}
$input_handle -> close;

#delete excess strings in sample name
my %hash_ss;
my @all_sample_name;
foreach my $every_id (keys %hash){
    my @one = @{$hash{$every_id}};
    push @all_sample_name,$one[0];
    $hash_ss{$every_id} = $one[0];
}

my %delete_result; #seqid => samplename
@all_sample_name = sort @all_sample_name;
for (my $i = 0; $i <= $#all_sample_name; $i += 2){
    my @one_word = split /$split_word/,$all_sample_name[$i];
    my @two_word = split /$split_word/,$all_sample_name[$i+1];
    die "@one_word\t@two_word\n" if scalar @one_word ne scalar @two_word;
    my @word = ();
    for (my $s = 0; $s <= $#one_word; $s ++){
        push @word,$one_word[$s] if $one_word[$s] eq $two_word[$s];
    }
    my $find_result = find_sequncingid(\%hash_ss,$all_sample_name[$i],$all_sample_name[$i+1]);
    my @find_result = @{$find_result};
    foreach my $every_result (@find_result){
        $delete_result{$every_result} = join("$split_word",@word);
    }
}

#create sample hash
my %sample;
foreach my $every (values %delete_result){
    my @library = grep { $delete_result{$_} eq $every } keys %delete_result; #length eq 2
    my @one_info = @{$hash{$library[0]}};
    my @two_info = @{$hash{$library[1]}};
    if ($one_info[1] =~ /linked/i){
        $sample{$every}{'UMILibraryR1'} = $one_info[2];
        $sample{$every}{'UMILibraryR2'} = $one_info[3];
        $sample{$every}{'AssembleLibraryR1'} = $two_info[2];
        $sample{$every}{'AssembleLibraryR2'} = $two_info[3];
    }else{
        $sample{$every}{'UMILibraryR1'} = $two_info[2];
        $sample{$every}{'UMILibraryR2'} = $two_info[3];
        $sample{$every}{'AssembleLibraryR1'} = $one_info[2];
        $sample{$every}{'AssembleLibraryR2'} = $one_info[3];
    }
}

print $output_handle "#SampleName\tUMILibraryR1\tUMILibraryR2\tAssembleLibraryR1\tAssembleLibraryR2\n";

foreach my $every (sort {$a cmp $b } keys %sample){
    print $output_handle "$every\t$sample{$every}{'UMILibraryR1'}\t$sample{$every}{'UMILibraryR2'}\t$sample{$every}{'AssembleLibraryR1'}\t$sample{$every}{'AssembleLibraryR2'}\n";
}
$output_handle -> close;

sub find_sequncingid {
    my %find_data = %{ shift @_ };
    my $one_w = shift;
    my $two_w = shift;
    my @result;
    #one
    foreach my $every_seqid (keys %find_data){
        push @result,$every_seqid if $find_data{$every_seqid} eq $one_w;
    }
    #two
    foreach my $every_seqid (keys %find_data){
        push @result,$every_seqid if $find_data{$every_seqid} eq $two_w;
    }
    die "$!" if scalar @result ne 2;
    return \@result;
}
