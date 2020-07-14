#!/usr/bin/perl
=head
    FileName: 05.TakeOut_UMI.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    05.TakeOut_UMI.pl : Part of the FASAS basic analysis script
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
    Input: Read-tag library sequence file
    Output: StepOne fastq
=cut

use strict;
use warnings;
use threads;
use POSIX qw(strftime);
use IO::File;
use apackages;
use File::Path qw(make_path remove_tree);
use Getopt::Long;
use File::Basename;

my ($fq,$rq,$type,$help,$rq_outputfile,$fq_outputfile);
my ($seq1,$seq2,$Regular,$Regular_rq,$r2_seq_length,$umi_length);
my ($alefta,$aleftp,$arighta,$arightp);
my $max_threads = 10;

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
    "f=s" => \$fq,
    "r=s" => \$rq,
    "t=s" => \$type,
    "aleftp=s" => \$aleftp,
    "alefta=s" => \$alefta,
    "arightp=s" => \$arightp,
    "arighta=s" => \$arighta,
    "umi_length=i" => \$umi_length,
    "max_threads=i" => \$max_threads,
    "help" => \$help,
) or echo_error("TakeOutUMI","Failed to get options");

if ($help or $option_number == 0){
    print "Take Out UMI Program
Usage:
    perl $0 -f [file] -r [file] -t [l|r]
Options:
    f -fq forword fastq[str,default:null]
    r -rq reverse fastq[str,default:null]
    t -type sequence type[str,l|r]
    h -help print this message
";
    exit;
}

my $type_word = "TakeOutUMI-".uc($type);

#check options
echo_error($type_word,"Forword Fastq: $fq not found") if ! -f $fq;
echo_error($type_word,"Reverse Fastq: $rq not found") if ! -f $rq;

chomp($alefta,$aleftp,$arighta,$arightp);
#get length
my $left_length = (length $aleftp) + $umi_length + (length $alefta);
my $right_length = (length $arightp) + $umi_length + (length $arighta);

#thread ref
my ($alefta_tr,$aleftp_tr,$arighta_tr,$arightp_tr);
$alefta_tr = reverse_complement($alefta);
$aleftp_tr = reverse_complement($aleftp);
$arighta_tr = reverse_complement($arighta);
$arightp_tr = reverse_complement($arightp);

$alefta_tr = ThreadRef($alefta_tr);
$aleftp_tr = ThreadRef($aleftp_tr);
$arighta_tr = ThreadRef($arighta_tr);
$arightp_tr = ThreadRef($arightp_tr);

$aleftp = ThreadRef($aleftp);
$alefta = ThreadRef($alefta);
$arightp = ThreadRef($arightp);
$arighta = ThreadRef($arighta);

my %ref_adapter = (
    'alefta' => $alefta,
    'alefta_tr' => $alefta_tr,
    'arighta' => $arighta,
    'arighta_tr' => $arighta_tr,
    'aleftp' => $aleftp,
    'aleftp_tr' => $aleftp_tr,
    'arightp' => $arightp,
    'arightp_tr' => $arightp_tr,
);

if ($type eq 'l'){
    $r2_seq_length = $left_length;
    $seq1 = $ref_adapter{'aleftp'};
    $seq2 = $ref_adapter{'alefta'};
    $rq_outputfile = 'read-tag_left_R2_one.fastq.gz';
    $fq_outputfile = 'read-tag_left_R1_one.fastq.gz';
    my @Regular = map { $_ =~ /right/ ? $ref_adapter{$_} : () } keys %ref_adapter;
    #added alefta because: error UMI seq
    $Regular = join('|',@Regular,$ref_adapter{'alefta'},$ref_adapter{'aleftp'});
    #added alefta because: error UMI seq
    $Regular_rq = join('|',$ref_adapter{'alefta'},$ref_adapter{'arightp_tr'},$ref_adapter{'aleftp_tr'},$ref_adapter{'arightp'},$ref_adapter{'alefta_tr'});
}elsif($type eq 'r'){
    $r2_seq_length = $right_length;
    $seq1 = $ref_adapter{'arightp'};
    $seq2 = $ref_adapter{'arighta'};
    $rq_outputfile = 'read-tag_right_R2_one.fastq.gz';
    $fq_outputfile = 'read-tag_right_R1_one.fastq.gz';
    my @Regular = map { $_ =~ /left/ ? $ref_adapter{$_} : () } keys %ref_adapter;
    #delete this seq if RIGHT_RI hava arighta or arightp
    $Regular = join('|',@Regular,$ref_adapter{'arighta'},$ref_adapter{'arightp'});
    #added arighta because: error UMI seq
    $Regular_rq = join('|',$ref_adapter{'arighta'},$ref_adapter{'arightp_tr'},$ref_adapter{'aleftp_tr'},$ref_adapter{'aleftp'},$ref_adapter{'arighta_tr'});
}else{
    echo_error($type_word,"Unkown type: $type");
}

#mkdir
my $local_time = strftime "%Y%m%d_%H:%M:%S", localtime;
my $tmpdir = "/dev/shm/get_UMI_${type}_$local_time";
make_path($tmpdir);

#readfile;
my $fq_handle = new IO::File "gzip -dc $fq |" or echo_error($type_word,"$!");
my $rq_handle = new IO::File "gzip -dc $rq |" or echo_error($type_word,"$!");

my @threads;
my $thread_number = 0;
my %data;
while (<$rq_handle>){
    chomp;
    my @header = split /\s+/,$_;
    my $seq = <$rq_handle>;
    my $middle = <$rq_handle>;
    my $qual = <$rq_handle>;
    my $header_fq = <$fq_handle>;
    my $seq_fq = <$fq_handle>;
    my $middle_fq = <$fq_handle>;
    my $qual_fq = <$fq_handle>;
    chomp($seq,$middle,$qual,$seq_fq,$middle_fq,$qual_fq);
    #
    next if length $seq < 100;
    $data{$header[0]} = "$seq\t$middle\t$qual\t$seq_fq\t$middle_fq\t$qual_fq";
    if (scalar keys %data == 500000){
        push @threads,threads -> create(\&ThreadSeq,\%data,$thread_number);
        $thread_number ++;
        %data = ();
        if (scalar @threads == $max_threads){
            my $every_thread = shift @threads;
            $every_thread -> join();
        }
    }
}
$fq_handle -> close;
$rq_handle -> close;

#analysis last data
if (%data){
    push @threads,threads -> create(\&ThreadSeq,\%data,$thread_number);
}

#join all threads
while (@threads){
    my $last_thread = shift @threads;
    $last_thread -> join();
}

my $result_file = `find $tmpdir -type f`;
chomp $result_file;
my @result_file = split "\r?\n",$result_file;
my @log_result_file = grep  { $_ =~ /log$/ } @result_file;
@log_result_file = sort { ($a =~ /(\d+)\.log$/)[0] <=> ($b =~ /(\d+)\.log$/)[0] } @log_result_file;
my $log_result_file = join(" ",@log_result_file);

my @fq_result_file = grep { $_ =~ /_R1_/ } @result_file;
@fq_result_file = sort { ($a =~ /(\d+)$/)[0] <=> ($b =~ /(\d+)$/)[0] } @fq_result_file;
my $fq_result_file = join(" ",@fq_result_file);

my @rq_result_file = grep { $_ =~ /_R2_/ } @result_file;
@rq_result_file = sort { ($a =~ /(\d+)$/)[0] <=> ($b =~ /(\d+)$/)[0] } @rq_result_file;
my $rq_result_file = join(" ",@rq_result_file);

my $last_log_file = basename($0,'.pl').'_'.$type.'.log';

system("cat $log_result_file > $last_log_file");
my $fq_strings = gzip_support($fq_outputfile,'output',4);
my $rq_strings = gzip_support($rq_outputfile,'output',4);
system("cat $fq_result_file $fq_strings");
system("cat $rq_result_file $rq_strings");

remove_tree($tmpdir);

#sub functions
sub CheckDupAd {
    my $sequence = shift;
    foreach my $every_adapter (values %ref_adapter){
        my $ex_number = () = $sequence =~ /$every_adapter/g;
        $ex_number >= 2 ? return 1 : next;
    }
}

sub ThreadSeq {
    my $th_fastq = shift;
    my $th_number = shift;
    my %th_fastq = %{$th_fastq};
    #output handle
    my $log_file = basename($0,'.pl').'_'.$type.'_'.$th_number.'.log';
    my $log_handle = new IO::File ">$tmpdir/$log_file" or die "$!";
    my $fq_outputfile_th = $fq_outputfile.'.'.$th_number;
    my $rq_outputfile_th = $rq_outputfile.'.'.$th_number;
    my $fq_output_handle = new IO::File ">$tmpdir/$fq_outputfile_th" or die "$!";
    my $rq_output_handle = new IO::File ">$tmpdir/$rq_outputfile_th" or die "$!";
    foreach my $line_fastq (keys %th_fastq){
        my @th_fastq = split "\t",$th_fastq{$line_fastq};
        if (CheckDupAd($th_fastq[0]) == 1 or CheckDupAd($th_fastq[3]) == 1){
            #check all sequence
            #if the sequence matches two or more common sequence
            #delete it
            print $log_handle "SEQDUP: $line_fastq\n$th_fastq[0]\n$th_fastq[3]\n";
            print $log_handle "---------------\n";
            next;
        }
        #check R1 sequence
        if ($th_fastq[3] =~ /$Regular/){
            print $log_handle "FQERROR: $line_fastq\n$th_fastq[3]\n";
            print $log_handle "---------------\n";
            next;
        }
        if ($th_fastq[0] =~ /^$seq1([ATCGN]{$umi_length})$seq2(.*)/){
            my $umi = $1;
            my $now_seq = $2;
            if ($now_seq =~ /$Regular_rq/){
                print $log_handle "RQERROR: $line_fastq\n$now_seq\n";
                print $log_handle "---------------\n";
                next;
            }
            my $now_qual = substr $th_fastq[2],$r2_seq_length; #N -> end
            print $fq_output_handle "$line_fastq\n$th_fastq[3]\n$th_fastq[4]\n$th_fastq[5]\n";
            print $rq_output_handle "$line_fastq $umi\n$now_seq\n$th_fastq[1]\n$now_qual\n";
        }elsif($th_fastq[0] =~ /^$seq1([ATCGN]{$umi_length})/){
            my $umi = $1;
            my $save_length = $r2_seq_length + 2;
            my $now_seq = substr $th_fastq[0],$save_length; #N -> end
            if ($now_seq =~ /$Regular_rq/){
                print $log_handle "RQERROR: $line_fastq\n$now_seq\n";
                print $log_handle "---------------\n";
                next;
            }
            my $now_qual = substr $th_fastq[2],$save_length; #N -> end
            print $fq_output_handle "$line_fastq\n$th_fastq[3]\n$th_fastq[4]\n$th_fastq[5]\n";
            print $rq_output_handle "$line_fastq $umi\n$now_seq\n$th_fastq[1]\n$now_qual\n";
        }
    }
    $fq_output_handle -> close;
    $rq_output_handle -> close;
    $log_handle -> close;
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
