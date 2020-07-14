#!/usr/bin/perl

=head
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Copyright 2018-2020 CapitalBio Corporation
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=cut

=function
    In this pm file, Sub-functions common to the PERL script are declared.
=cut

package apackages;
require Exporter;

use strict;
use warnings;
use Switch;
use IO::File;
use Cwd qw(abs_path);
use File::Which qw(which);
use POSIX qw(strftime);

our @ISA = qw(Exporter);
our @EXPORT = qw(echo_error echo_info echo_warn CheckCPU CheckMem reverse_complement gzip_support read_complex_fast guess_fast_format);

sub CheckCPU {
    #1. monitor CPU usage for one second
    #2. calculate available CPU resources and add one
    #3. returns the number of available CPU threads
    my $cpuinfo_handle = new IO::File "/proc/cpuinfo" or die "$!";
    my @cpuinfo = <$cpuinfo_handle>;
    @cpuinfo = grep { $_ =~ /^processor\s+:\s+\d+/ } @cpuinfo;
    my $cpu_number = scalar @cpuinfo;
    $cpuinfo_handle -> close;
    my $stat_handle = new IO::File "/proc/stat" or die "$!";
    my $stat_one_line = <$stat_handle>;
    chomp $stat_handle;
    my @stat_one_line = split /\s+/,$stat_one_line;
    $stat_handle -> close;
    sleep 1;
    $stat_handle = new IO::File "/proc/stat" or die "$!";
    $stat_one_line = <$stat_handle>;
    chomp $stat_handle;
    my @stat_two_line = split /\s+/,$stat_one_line;
    my $jiff_zero = 0;
    my $jiff_one = 0;
    for (my $i = 1; $i <= 7; $i ++){
        $jiff_zero += $stat_one_line[$i];
        $jiff_one += $stat_two_line[$i];
    }
    my $sys_idle = ($stat_one_line[4] - $stat_two_line[4]) / ($jiff_zero - $jiff_one);
    my $sys_idle_number = int($cpu_number * $sys_idle);
    return $sys_idle_number;
}

sub CheckMem {
    #MemTotal = MemFree +[Slab+ VmallocUsed + PageTables + KernelStack + HardwareCorrupted + Bounce + X]+[Cached + AnonPages + Buffers + (HugePages_Total * Hugepagesize)])
    my $meminfo_file = '/proc/meminfo';
    if ( -e $meminfo_file ){
        my $available_mem_size = 0; #units is GB
        my $handle_meminfo = new IO::File "$meminfo_file" or die "$!";
        local $/ = undef;
        my $meminfo = <$handle_meminfo>;
        if ($meminfo =~ /MemAvailable:\s+(\d+)\s+(\w+)/){
            my $available_mem_value = $1;
            my $available_mem_units = $2;
            switch ($available_mem_units) {
                case /kb/i { $available_mem_size = $available_mem_value/1024/1024 }
                case /mb/i { $available_mem_size = $available_mem_value/1024 }
                case /gb/i { $available_mem_size = $available_mem_value }
                else { return -2 }
            }
            $handle_meminfo -> close;
            return $available_mem_size;
        }else{
            my $available_mem_size = 0;
            my $MemFree = 0;
            my $Buffer = 0;
            my $Cached = 0;
            my $units = '';
            if ($meminfo =~ /MemFree:\s+(\d+)\s+(\w+)/){
                $MemFree = $1;
                $units = $2;
            }
            if ($meminfo =~ /Buffer:\s+(\d+)\s+\w+/){
                $Buffer = $1;
            }
            if ($meminfo =~ /Cached:\s+(\d+)\s+\w+/){
                $Cached = $1;
            }
            my $available_mem_value = $MemFree + $Buffer + $Cached;
            switch ($units) {
                case /kb/i { $available_mem_size = $available_mem_value/1024/1024 }
                case /mb/i { $available_mem_size = $available_mem_value/1024 }
                case /gb/i { $available_mem_size = $available_mem_value }
                else { return -2 }
            }
            $handle_meminfo -> close;
            return $available_mem_size;
        }
    }else{
        return -1;
    }
}

sub echo_info {
    my $local_time = strftime "%Y%m%d-%H:%M:%S", localtime;
    print "[$local_time][$_[0]]: <Info> $_[1]\n";
}

sub echo_error {
    my $local_time = strftime "%Y%m%d-%H:%M:%S", localtime;
    print "[$local_time][$_[0]]: <Error> $_[1]\n";
    exit 1;
}

sub echo_warn {
    my $local_time = strftime "%Y%m%d-%H:%M:%S", localtime;
    print "[$local_time][$_[0]]: <Warn> $_[1]\n";
}

sub reverse_complement {
    my $input = shift;
    $input =~ tr/ATCGRYMKSWHBVDN/TAGCYRKMWSDVBHN/;
    $input = reverse $input;
    return $input;
}

sub gzip_support {
    my $file = shift; #file address
    my $type = shift; #input or output or append output
    my $max_threads = shift;
    if ($file =~ /gz$/ or $file =~ /gzip$/){
        my $use_command = '';
        if (which('pigz')){
            $use_command = 'pigz';
        }elsif (which('gzip')){
            $use_command = 'gzip';
        }else{
            echo_error("gzip_support","gzip and pigz not found!");
        }
        #max_threads
        if ($use_command eq 'pigz' and defined $max_threads){
            $use_command .= ' -p '.$max_threads;
        }
        my $key_strings = '';
        switch ($type) {
            case 'output' { $key_strings = "| $use_command > $file"; }
            case 'append' { $key_strings = "| $use_command >> $file"; }
            case 'input' { $key_strings = "$use_command -dc $file |"; }
            else { echo_error("gzip_support","unknown type: $type !"); }
        }
        return $key_strings;
    }else{
        my $key_strings = '';
        switch ($type) {
            case 'output' { $key_strings = "> $file"; }
            case 'append' { $key_strings = ">> $file"; }
            case 'input' { $key_strings = "$file"; }
            else { echo_error("gzip_support","unknown type: $type !"); }
        }
        return $key_strings;
    }
}

sub trap_function {
    #options: script_name and temp_folder
    #When this function is triggered, it will clean up the temporary directory to ensure system resources are available.
    my $script_name = shift;
    my $dir = shift;
    echo_info($script_name,"This script receives the termination signal and is exiting");
    my @dir = @{ $dir };
    foreach my $every_dir (@dir) {
        my $abs_folder_path = abs_path($every_dir);
        if ( -e $abs_folder_path ){
            unlink $abs_folder_path;
        }else{
            my $warn_words = 'File or folder not found: '.$abs_folder_path;
            echo_warn($script_name,$warn_words);
        }
    }
    exit;
}

#annotation
our $read_fast_register = '#orgin#';
#
sub read_complex_fast {
    my $input_handle = shift;
    my ($id,$decs);
    my $seq = '';
    if ($read_fast_register ne '#orgin#'){
        ($id, $decs) = _thread_fast_head($read_fast_register);
    }else{
        my $oneline = <$input_handle>;
        chomp $oneline;
        ($id, $decs) = _thread_fast_head($oneline);
    }
    while (my $line = <$input_handle>){
        chomp $line;
        $read_fast_register = $line and last if $line =~ /^[>@]/;
        $seq .= $line;
    }
    #return undef if eof $input_handle;
    return $id,$decs,$seq;
}

sub _thread_fast_head {
    my $head_strings = shift;
    if ($head_strings =~ /^[>@](\S+?)\s(.+)$/){
        return $1,$2;
    }elsif ($head_strings =~ /^[>@](\S+?)$/){
        return $1,'none';
    }else{
        return 'none','none';
    }
}

sub guess_fast_format {
    my $input_file = shift;
    my @data;
    my $input_handle = new IO::File "$input_file" or echo_error("Guess_fast_format","$!");
    @data[0..4] = <$input_handle>;
    if ($data[0] =~ /^>/){
        if ($data[2] =~ /^>/){
            return 'simple_fasta';
        }elsif($data[2] =~ /^[a-zA-Z]+$/){
            return 'complex_fasta';
        }else{
            return 'bad_fasta';
        }
    }elsif ($data[0] =~ /^@/){
        if ($data[4] =~ /^@/){
            return 'simple_fastq';
        }elsif($data[4] =~ /^[a-zA-Z]+$/){
            return 'complex_fastq';
        }else{
            return 'bad_fastq';
        }
    }else{
        return 'bad_format';
    }
}
