#!/usr/bin/perl
=head
    FileName: Create_CofigFile.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Create_ConfigFile.pl
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
    Input: Sample information table
    Output: FASAS ConfigFile for every sample
=cut

=input_data
#Comment line
#SampleName\tUMILibraryR1\tAssembleLibraryR1
or
#SampleName\tUMILibraryR1\tUMILibraryR2\tAssembleLibraryR1\tAssembleLibraryR2
=cut

=modify_hash
default options are in the 'default_parameter' hash
change it before using this script
=cut

use strict;
use warnings;
use IO::File;
use Switch;
use Getopt::Long;

my ($input_tsv,$output_dir,$help);
my %default_parameter = (
    'SampleName' => 'Test',
    'WorkFolder' => './',
    'LibraryInfo' => '/inspurfs/transbio_sx02/ZhangKe/codes/16sFASAP/data/LibraryInfo/126_length_library_info.txt',
    'UMILibraryR1' => '',
    'UMILibraryR2' => '',
    'AssembleLibraryR1' => '',
    'AssembleLibraryR2' => '',
    'CoverageDatabase' => '/inspurfs/database/metagenome/16SMicrobial_uniq/16SMicrobial',
    'PrimerDatabaseFasta' => '/inspurfs/database/metagenome/16SMicrobial_uniq/16SMicrobial',
    'PrimerDatabaseTaxonomy' => '/inspurfs/database/metagenome/16SMicrobial_uniq/16SMicrobial_taxonomy.tsv',
    'AssembleProgram' => 'idba_ud',
    'AssembleIlluminaAdapter' => 'CTGTCTCTTATACACATCT',
    'ContigLength' => 1200,
);

my $options_number = scalar @ARGV;
GetOptions (
    "input_tsv=s" => \$input_tsv,
    "output_dir=s" => \$output_dir,
    "help" => \$help,
) or die "$!";

if ($help or $options_number == 0){
    print "Config File Create Program
usage:
    !!modify the 'default_parameter' hash in the script before using this script!!
    perl $0 --input_tsv [file] --output_dir [path]
options:
    input_tsv     strings    a input file, including information about the sample name and sequence file location [required]
    output_dir    strings    a folder, storing Config File [required]
    help          switch     get this information [off]
";
    exit;
}

my @parameter_array = qw(SampleName WorkFolder LibraryInfo UMILibraryR1 UMILibraryR2 AssembleLibraryR1 AssembleLibraryR2 CoverageDatabase PrimerDatabaseFasta PrimerDatabaseTaxonomy AssembleProgram AssembleIlluminaAdapter ContigLength);

my $input_handle = new IO::File "$input_tsv" or die "$!";
mkdir $output_dir if ! -e $output_dir;

my $line_index = 0;
while (<$input_handle>){
    chomp;
    $line_index ++;
    next if $_ =~ /^#/ or $_ =~ /^$/;
    my @line = split "\t",$_;
    switch ($#line){
        case 2 { ThreadLine(\@line,3,$output_dir) }
        case 4 { ThreadLine(\@line,5,$output_dir) }
        else { WarnLine{$#line,$line_index} }
    }
}
close $input_handle;

sub ThreadLine {
    my @line = @{ shift @_ };
    my $col_number = shift;
    my $output_path = shift;
    my $config_file = $output_path.'/'.$line[0].".ConfigFile";
    my $output_handle = new IO::File ">$config_file" or die "$!";
    if ($col_number == 3){
        foreach my $every_array (@parameter_array){
            switch ($every_array) {
                case 'SampleName' {print $output_handle "$every_array = $line[0]\n";}
                case 'WorkFolder' {print $output_handle "$every_array = ./$line[0]\n";}
                case 'UMILibraryR1' {print $output_handle "$every_array = $line[1]\n";}
                case 'UMILibraryR2' {PrintAdd($line[1],'UMILibraryR2',$output_handle)}
                case 'AssembleLibraryR1' {print $output_handle "$every_array = $line[2]\n";}
                case 'AssemlbeLibraryR2' {PrintAdd($line[2],'AssemlbeLibraryR2',$output_handle)}
                else {print $output_handle "$every_array = $default_parameter{$every_array}\n";}
            }
        }
    }elsif ($col_number == 5 ){
        foreach my $every_array (@parameter_array){
            switch ($every_array) {
                case 'SampleName' {print $output_handle "$every_array = $line[0]\n";}
                case 'WorkFolder' {print $output_handle "$every_array = ./$line[0]\n";}
                case 'UMILibraryR1' {print $output_handle "$every_array = $line[1]\n";}
                case 'UMILibraryR2' {print $output_handle "$every_array = $line[2]\n";}
                case 'AssembleLibraryR1' {print $output_handle "$every_array = $line[3]\n";}
                case 'AssembleLibraryR2' {print $output_handle "$every_array = $line[4]\n";}
                else {print $output_handle "$every_array = $default_parameter{$every_array}\n";}
            }
        }
    }
    $output_handle -> close;
}

sub WarnLine {
    my $one_line_number = shift;
    my $line_index = shift;
    $one_line_number ++;
    die "CreateConfig Error: The number of columns in the ${line_index}th row is $one_line_number, it looks wrong.\n";
}

sub PrintAdd {
    my $R1_path = shift;
    my $type = shift;
    my $output_handle = shift;
    $R1_path =~ s/(.*)R1(.+?)$/$1R2$2/;
    print $output_handle "$type = $R1_path\n";
}
