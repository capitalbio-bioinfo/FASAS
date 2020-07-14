#!/usr/bin/perl
=head
    FileName: 08.Parallel_Assemble.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    08.Parallel_Assemble.pl : Part of the FASAS basic analysis script
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
    Input: All 16s rRNA sequence
    Output: 16s rRNA cnotig
=cut

use warnings;
use strict;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use Switch;
use IO::File;
use apackages;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Parallel::ForkManager;

#default options
my ($input_file,$input_file2,$assemble_program,$help);
my $max_thread = 40;
my $output_file = './contigs.fa';

#get assemble script path
my $script_path = abs_path($0);
my $script_folder_path = dirname($script_path);
my $assemble_script = $script_folder_path.'/../bin/run_assemble.sh';
echo_error("ParallelAssemble","No assemble script, $assemble_script") if ! -e $assemble_script;

my $option_number = scalar @ARGV;
GetOptions (
    "forward_fastq=s" => \$input_file,
    "reverse_fastq=s" => \$input_file2,
    "output_file=s" => \$output_file,
    "max_thread=s" => \$max_thread,
    "assemble_program=s" => \$assemble_program,
    "help" => \$help,
) or echo_error("ParallelAssemble","Failed to get options");

#help
if ($help or $option_number <= 0){
    print "Parallel Assemble of Full-Length 16s data Program
Usage:
    forward_fastq    path     Input file, forward FASTQ file [required]
    reverse_fastq    path     Input file, reverse FASTQ file [required]
    output_file      strings  Output file, assemble result, FASTA format [required]
    max_thread       int      Maximum number of threads [40]
    assemble_program strings  Choose a program to assemble [cap3 or idba_ud]
    help             switch   Get this help information [off]

";
    exit;
}

#check options
echo_error("ParallelAssemble","$input_file cannot be read!") if ! -r $input_file;
echo_error("ParallelAssemble","$input_file2 cannot be read!") if ! -r $input_file2;
echo_error("ParallelAssemble","No support assembly program $assemble_program") if $assemble_program !~ /cap3/i and $assemble_program !~ /idba_ud/i;
$assemble_program = lc($assemble_program);
echo_error("ParallelAssemble","Failed to find assemble script!") if ! -r $assemble_script;

#open handle
my $input_r1_strings = gzip_support($input_file,'input',2);
my $input_r2_strings = gzip_support($input_file2,'input',2);
my $input_fastq_handle = new IO::File "$input_r1_strings" or echo_error("ParallelAssemble","$!");
my $input_fastq_handle2 = new IO::File "$input_r2_strings" or echo_error("ParallelAssemble","$!");
my $log_handle = new IO::File ">parallel_assemble_dbug.log" or echo_error("ParallelAssemble","$!");

#read fastq
echo_info("ParallelAssemble","Read FASTQ file");
my %nameseq;
my %share_key;
while (<$input_fastq_handle>){
    chomp;
    my $seq = <$input_fastq_handle>;
    my $mid = <$input_fastq_handle>;
    my $qual = <$input_fastq_handle>;
    my $header2 = <$input_fastq_handle2>;
    my $seq2 = <$input_fastq_handle2>;
    my $mid2 = <$input_fastq_handle2>;
    my $qual2 = <$input_fastq_handle2>;
    if ($_ =~ /@(.+)_\d+/){
        $nameseq{$1} .= ">${_}\n".$seq.">$header2".$seq2;
        $share_key{$1} = 1;
    }else{
        echo_warn("ParallelAssemble","Regular expression match failed while reading in sequence\n$_");
    }
}
$input_fastq_handle -> close;
$input_fastq_handle2 -> close;

echo_info("ParallelAssemble","Create TEMP folder");
#better IOPS with shm
sleep 1; #Avoid duplicate folder names
my $temp_dirname;
my $local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
my $tmp_name_before = '';
$assemble_program =~ /idba_ud/i ? ( $tmp_name_before = 'IDBA_' ) : ( $tmp_name_before = 'Cap3_' );
if (-w "/dev/shm"){
    $temp_dirname = '/dev/shm/'.$tmp_name_before.$local_time;
    echo_info("ParallelAssemble","using shm");
}else{
    $temp_dirname = '/tmp/'.$tmp_name_before.$local_time;
    echo_info("ParallelAssemble","using hdd");
}

mkdir $temp_dirname or echo_error("ParallelAssemble","Failed to generate TEMP folder");

#parallel ready
echo_info("ParallelAssemble","Parallel Ready");
foreach my $every_key (keys %nameseq){
    my $p_dir = $temp_dirname.'/'.$every_key;
    mkdir $p_dir or die "$!";
    open OUT,">","$p_dir/$every_key.fa" or die "$!";
    print OUT $nameseq{$every_key};
    close OUT;
}
my $available_cpu_number = CheckCPU() + 1;
$max_thread = $available_cpu_number if $available_cpu_number < $max_thread;

#clean nameseq
undef %nameseq;

echo_info("ParallelAssemble","Parallel Assembly");
my $pm = Parallel::ForkManager->new($max_thread);
BACKEND:
foreach my $every_key (keys %share_key){
    $pm->start and next BACKEND;
    my $p_dir = $temp_dirname.'/'.$every_key;
    my $stat = system("bash $assemble_script --fasta_file $p_dir/$every_key.fa --output_folder $p_dir --program $assemble_program");
    if ($stat ne 0){
        print $log_handle "$every_key: Assembe Script return $stat\n";
    }
    $pm -> finish;
}
$pm -> wait_all_children;

echo_info("ParallelAssemble","Read Contigs");
my $umi_dir = `find $temp_dirname -mindepth 1 -maxdepth 1 -type d`;
my @umi_dir = split "\r?\n",$umi_dir;

my $handle_out = new IO::File ">$output_file" or echo_error("ParallelAssemble","$output_file $!");
foreach my $every_dir (@umi_dir){
    switch ($assemble_program) {
        case /idba_ud/i { ReadIDBAudResult($every_dir) }
        case /cap3/i { ReadCAP3Result($every_dir) }
        else { echo_error("ParallelAssemble","Unknown assembly program $assemble_program in $every_dir"); }
    }
}
$handle_out -> close;
$log_handle -> close;

#delete tmp folder
if ($assemble_program =~ /cap3/){
    remove_tree($temp_dirname) or echo_error("ParallelAssemble","Failed to Remove Cap3 Folder");
}else{
    my $temp_handle = new IO::File ">./idba_ud_address.txt" or echo_error("ParallelAssemble","$!");
    print $temp_handle "$temp_dirname";
    $temp_handle -> close;
}

sub ReadCAP3Result {
    my $every_dir = shift;
    my $max_length = 0;
    my $long_contig = '';
    my $file_basename = basename($every_dir);
    my $dirname = $file_basename.'_cap3_contigs';
    my $cap3_contig_file = $every_dir.'/'.$file_basename.'_r1.fa.cap.contigs';
    my $cap3_handle = new IO::File "$cap3_contig_file" or echo_error("ParallelAssemble","$!");
    while (my ($id, $decs, $seq) = read_complex_fast($cap3_handle)){
        if (length $seq >= $max_length){
            $max_length = length $seq;
            $long_contig = $seq;
        }
        last if eof($cap3_handle);
    }
    print $handle_out ">$dirname\n$long_contig\n";
}

sub ReadIDBAudResult {
    my $every_dir = shift;
    my $kmer_ref = ReadIDBALog($every_dir);
    echo_warn("ParallelAssemble","In $every_dir, IDBA_ud LogFile not found") if $kmer_ref =~ /logfile not found/;
    my %kmer = %{$kmer_ref};
    my $good_kmer = 0;
    my $good_contig_number = 0;
    my $good_max_contig_length = 0;
    my $good_total_length = 0;
    #find best kmer
    foreach my $every_kmer (sort {$a <=> $b} keys %kmer){
        #print "middle: $good_kmer\t$good_contig_number\t$good_max_contig_length\t$good_total_length\n";
        if ($good_kmer == 0){
            $good_kmer = $every_kmer;
            $good_contig_number = @{$kmer{$every_kmer}}[0];
            $good_max_contig_length = @{$kmer{$every_kmer}}[1];
            $good_total_length = @{$kmer{$every_kmer}}[2];
        }else{
            #new max contig length > old max contig length
            if (@{$kmer{$every_kmer}}[1] > $good_max_contig_length){
                $good_kmer = $every_kmer;
                $good_contig_number = @{$kmer{$every_kmer}}[0];
                $good_max_contig_length = @{$kmer{$every_kmer}}[1];
                $good_total_length = @{$kmer{$every_kmer}}[2];
            }
            #new max contig length >= old max contig length && new contig number < old contig number
            if (@{$kmer{$every_kmer}}[1] >= $good_max_contig_length and @{$kmer{$every_kmer}}[0] < $good_contig_number){
                $good_kmer = $every_kmer;
                $good_contig_number = @{$kmer{$every_kmer}}[0];
                $good_max_contig_length = @{$kmer{$every_kmer}}[1];
                $good_total_length = @{$kmer{$every_kmer}}[2];
            }
        }
    }
    #print "result: $good_kmer\t$good_contig_number\t$good_max_contig_length\t$good_total_length\n";
    #read idba_ud contig
    my $contigs_file = $every_dir.'/contig-'.$good_kmer.'.fa';
    my $max_length = 0;
    my $long_contig = '';
    my $dirname = basename(dirname($contigs_file));
    $dirname .= "_${good_kmer}_${good_contig_number}_${good_max_contig_length}_${good_total_length}";
    my $idba_handle = new IO::File "$contigs_file" or echo_error("ParallelAssemble","$!");
    while (my ($id, $decs, $seq) = read_complex_fast($idba_handle)){
        if (length $seq >= $max_length){
            $max_length = length $seq;
            $long_contig = $seq;
        }
        last if eof($idba_handle);
    }
    print $handle_out ">$dirname\n$long_contig\n";
}

sub ReadIDBALog {
    my $logdir = shift;
    my %kmer_read;
    my $kmer;
    return 'logfile not found' if ! -e "$logdir/log";
    my $log_handle = new IO::File "<$logdir/log" or echo_error("ParallelAssemble","$logdir/log no read");
    while (<$log_handle>){
        chomp;
        if ($_ =~ /^kmer\s+(\d+)/){
            $kmer = $1;
        }elsif ($_ =~ /^contigs:/ and defined $kmer){
            my @contigs_info = split " ",$_;
            if ($#contigs_info < 10){
                echo_warn("$logdir, $kmer, contigs_info length is too short");
                undef $kmer;
                next;
            }
            $kmer_read{$kmer} = [$contigs_info[1],$contigs_info[5],$contigs_info[10]];
            undef $kmer;
        }
    }
    $log_handle -> close;
    return \%kmer_read;
}
