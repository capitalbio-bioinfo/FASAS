#!/usr/bin/perl
=head
    FileName: ktImportText.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    ktImportText.pl
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
    Output: Krona result
=cut

use strict;
use warnings;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use File::Basename;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use File::Which;
use File::Find;

my ($input_dir,$output_file,$help);

my $options_number = scalar @ARGV;
GetOptions (
    'input_dir=s' => \$input_dir,
    'output_file=s' => \$output_file,
    'help' => \$help,
);

if ($help or $options_number < 4){
    print "ktImportText Workflow Program
Usage:
    perl $0 -i [FL16S Project folder] -o [ktImportText outfile]
Options:
    input_dir      path       A full-length 16S project folder [required]
    output_file    strings    Output file of ktImportText [required]
    help           switch     Get this message [off]
";
    exit;
}

#check options
die "$input_dir not found!\n" if ! -d $input_dir;

#check ktImportText
my $local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
my $kyImportText_path = which 'ktImportText';
$kyImportText_path ? print "ktImportText_[$local_time]: Confirm that the ktImportText program exists!\n" : die "ktImportText_ERROR: ktImportText program not found!\n";

#search file
my @file;
$input_dir = abs_path($input_dir);
find(\&FindFile,$input_dir);

#tmp folder
$local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
my $tmpdir = '/tmp/ktImportText_'.$local_time.'_TmpDir';
make_path($tmpdir) ? print "KtImportText_[$local_time]: Create TMP Dir $tmpdir\n" : die "KtImportText_ERROR: Failed to Create TMP Dir \n";

#thread file
foreach my $every_file (@file){
    my $basename = basename($every_file,'_species.tsv');
    my $every_tmp_file = $tmpdir.'/'.$basename.'.tsv';
    my $every_file_handle = new IO::File "$every_file" or die "$!";
    my $every_tmp_handle = new IO::File ">$every_tmp_file" or die "$!";
    while (<$every_file_handle>){
        chomp;
        next if $_ =~ /^Tax/;
        my @line = split "\t",$_;
        $line[0] =~ s/;/\t/g;
        print $every_tmp_handle "$line[1]\t$line[0]\n";
    }
    $every_file_handle -> close;
    $every_tmp_handle -> close;
}

my $tmp_file = `find $tmpdir -type f |grep tsv`;
my @tmp_file = split "\r?\n",$tmp_file;
@tmp_file = sort @tmp_file;
my $krona_input = join(" ",@tmp_file);

my $krona_stat = system("ktImportText $krona_input -o $output_file");
$local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
$krona_stat == 0 ? print "KtImportText_[$local_time]: ktImportText Success\n" : die "KtImportText_ERROR: ktImportText return $krona_stat!\n";

remove_tree($tmpdir) ? print "Remove TMP Dir $tmpdir\n" :  die "KtImportText: Remove TMP Dir ERROR\n";

sub FindFile {
    if ( -f $File::Find::name ){
        if ( $File::Find::name =~ /species\.tsv$/ ){
            push @file,$File::Find::name;
        }
    }
}

sub FindTMP {
    if ( -f $File::Find::name ){
        if ( $File::Find::name =~ /tsv$/){
            push @tmp_file,$File::Find::name;
        }
    }
}
