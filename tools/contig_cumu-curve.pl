#!/usr/bin/perl
=head
    FileName: Contig_Len2Num.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Contig_Len2Num.pl
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
    Output: Cumulative curve of Contig length
=cut

use strict;
use warnings;
use apackages;
use Getopt::Long;
use IO::File;
use File::Find;
use File::Basename;
use Cwd qw(abs_path realpath);
use List::Util qw(uniq);

my ($input_folder,$output_file,$pdf_output,$reverse_filter,$sample_filter,$help);

#default opetions
our $file_name = 'fix_contigs.fa';

my $options_number = scalar @ARGV;
GetOptions (
    "input_folder=s" => \$input_folder,
    "output_file=s" => \$output_file,
    "sample_filter=s" => \$sample_filter,
    "reverse_filter" => \$reverse_filter,
    "pdf_output=s" => \$pdf_output,
    "file_name=s" => \$file_name,
    "help" => \$help,
) or echo_error("Contig_Len2Num","Failed to get options!");

if ($help or $options_number < 1){
    print "Cumulative Frequency of Contig Length Plotting Program
Usage:
    perl $0 --input_folder|-i [FASAS project] --output_file|-o [output] --pdf_output|-p [pdf]
Options:
    input_folder    strings    There are N FASAS analysis results in this directory [require].
    output_file     strings    The intermediate result of the program outputting the cumulative value during the calculation process [require].
    sample_filter   strings    filter sample, by regular expression [null].
    reverse_filter  switch     reverse filter sample [off].
    pdf_output      strings    Cumulative curve image of the final output of the program [require].
    file_name       strings    Contig file name [fix_contigs.fa].
    help            switch     print this messages [off].
";
    exit;
}

echo_error("Contig_Len2Num","Input folder not found!") if ! -d $input_folder;

our @file = ();
find(\&find_condition, $input_folder);
echo_error("Contig_Len2Num","No files were found!") if scalar @file < 1;

#sort out;
@file = map { abs_path($_) } @file;
@file = map { realpath($_) } @file;
@file = sort(uniq(@file));
if ($sample_filter){
    if ($reverse_filter){
        @file = grep { $_ !~ /$sample_filter/ } @file;
    }else{
        @file = grep { $_ =~ /$sample_filter/ } @file;
    }
}

my %stat;
foreach my $every_file (@file){
    my $project_dir = dirname(dirname($every_file));
    my $porject_name = basename($project_dir);
    my $input_handle = new IO::File "$every_file" or echo_error("Contig_Len2Num","Failed to open handle");
    while (<$input_handle>){
        my $seq = <$input_handle>;
        chomp $seq;
        my $length = length $seq;
        $stat{$porject_name}{$length} += 1;
    }
    $input_handle -> close;
}

my %cumu;
foreach my $every_project (sort {$a cmp $b} keys %stat){
    my $before_length = -1;
    foreach my $every_length (sort {$a <=> $b} keys %{$stat{$every_project}}){
        if ($before_length == -1){
            $cumu{$every_project}{$every_length} = $stat{$every_project}{$every_length};
            $before_length = $every_length;
        }else{
            $cumu{$every_project}{$every_length} = $stat{$every_project}{$every_length} + $cumu{$every_project}{$before_length};
            $before_length = $every_length;
        }
    }
}

my $output_handle = new IO::File ">$output_file" or echo_error("Contig_Len2Num","Failed to open handle");
print $output_handle "Sample\tLength\tFreq.\n";
foreach my $every_project (sort {$b cmp $a} keys %cumu){
    foreach my $every_length (sort {$a <=> $b} keys %{$cumu{$every_project}}){
        print $output_handle "$every_project\t$every_length\t$cumu{$every_project}{$every_length}\n";
    }
}
$output_handle -> close;

#plot
#Rscipt ${FASASHome}/r/tools-contig_cumu-curve.r $output_file $pdf_output
echoerr("Contig_Len2Num","FASASHome ENV not found!") if ! $ENV{'FASASHome'};
if (-e "$ENV{'FASASHome'}/r/tools-contig_cumu-curve.r"){
    my $rscript_stat = system("Rscript $ENV{'FASASHome'}/r/tools-contig_cumu-curve.r $output_file $pdf_output");
    if ($rscript_stat != 0){
        echo_error("Contig_Len2Num","Failed to Run R script!");
    }else{
        echo_info("Contig_Len2Num","Run R script sccuessful!");
    }
}else{
    echo_error("Contig_Len2Num","tools-contig_cumu-curve.r not found!");
}

sub find_condition {
    if ( $File::Find::name =~ /$file_name$/ ){
        push @file,$File::Find::name;
    }
}
