#!/usr/bin/env perl
=head
    FileName: Create_QsubScript.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Create_QsubScript.pl
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
    Input: FASAS ConfigFile folder
    Output: PBS script for every ConfigFile
=cut

=pbs
    #!/bin/bash
    #PBS -N [sample name]
    #PBS -l nodes=1:ppn=20
    #PBS -o [log_path]
    #PBS -e [err_path]
    date
    cd [run_path]
    bash [FASAS_path]run_analysis.sh -f [config_file]
    date
=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my ($log_path,$run_path,$config_path,$FASAS_name,$help);
$FASAS_name = 'run_analysis.sh';

my $options_number = scalar @ARGV;
GetOptions (
    'pbs_folder=s' => \$log_path,
    'run_path=s' => \$run_path,
    'config_folder=s' => \$config_path,
    'FASAS_name=s' => \$FASAS_name,
    'help' => \$help,
) or die "$!";

if ($help or $options_number < 1){
    print "Qsub Script Create Script
Usage:
    perl $0 -p|--pbs_folder [path] -r|--run_path [path] -c|--config_folder [path]
Options:
    pbs_folder       strings    Output, folder for pbs scripts [required]
    run_path         strings    Working folder of the pbs script [required]
    config_folder    strings    Folder of ConfigFile [required]
    FASAS_name       strings    FASAS main program name [run_analysis.sh]
    help             switch     Show this help message and exit [off]
";
    exit;
}

die "need pbs_folder\n" if ! defined $log_path;
die "need run_path\n" if ! defined $run_path;
die "need config_folder\n" if ! defined $config_path;

#find run_analysis
my $dirname = dirname(__FILE__);
my $run_analysis = $dirname.'../'.$FASAS_name;
$run_analysis = abs_path($run_analysis);

#get abs_path
$log_path = abs_path($log_path);
$log_path .= '/' if $log_path !~ /\/$/;
$run_path = abs_path($run_path);
$run_path .= '/' if $run_path !~ /\/$/;

my $config = `find $config_path -type f`;
my @config = split "\r?\n",$config;

#filter config file
@config = grep { $_ =~ /ConfigFile$/ }  @config;

die "No files were found" if scalar @config == 0;
foreach my $every (@config){
    my $every = abs_path($every);
    my $bname = basename($every,'.ConfigFile');
    my $e = $log_path.$bname.'.err';
    my $l = $log_path.$bname.'.log';
    my $pbs_file = $log_path.$bname.'.pbs';
    my $output_handle = new IO::File ">$pbs_file" or die "$!";
    print $output_handle "#!/bin/bash\n
#PBS -N $bname
#PBS -l nodes=1:ppn=20
#PBS -o $l
#PBS -e $e

date
cd $run_path
bash $run_analysis -f $every
date\n";
    $output_handle -> close;
}
