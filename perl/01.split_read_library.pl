#!/usr/bin/perl
=head
    FileName: 01.Split_Assemble_Library.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    01.Split_Assemble_Library.pl : Part of the FASAS basic analysis script
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
    Input: Linked-Tag library sequence
    Output: Aligned sequence
=cut

=date_structure
FASTQ gzip or no gzip
R2:
    @.+
    [regular sequence][ATCGN]+
    \+
    FFFFF+
=cut

use strict;
use warnings;
use IO::File;
use apackages;
use Getopt::Long;
use List::Util qw(min uniq);

my ($r1_file,$r2_file,$log_file,$alp,$arp,$max_distance,$help);

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

my %dege_out = (
    'R' => 'A',
    'Y' => 'C',
    'M' => 'A',
    'K' => 'G',
    'S' => 'G',
    'W' => 'A',
    'H' => 'A',
    'B' => 'G',
    'V' => 'G',
    'D' => 'G',
    'N' => 'A',
);

#default options
$log_file = 'read-tag_library.log';
$max_distance = 3;

my $options_number = scalar @ARGV;
GetOptions (
    "f=s" => \$r1_file,
    "r=s" => \$r2_file,
    "alp=s" => \$alp,
    "arp=s" => \$arp,
    "max_distance=i" => \$max_distance,
) or echo_error("SplitAssembleLibrary","Failed to get options");

if ($help or $options_number < 5){
    print "

";
    exit 0;
}

echo_error("SplitAssembleLibrary","$r1_file no exists!") if ! -e $r1_file;
echo_error("SplitAssembleLibrary","$r2_file no exists!") if ! -e $r2_file;

#thread ref sequence
##left
chomp $alp;
my $alp_length = length $alp;
my @alp = split "",$alp;
my @alp_out;
for (my $i = 0; $i <= $#alp; $i ++){
    if (exists $dege{$alp[$i]}){
        push @alp_out,$dege_out{$alp[$i]};
        $alp[$i] = $dege{$alp[$i]};
    }else{
        push @alp_out,$alp[$i];
    }
}
$alp = join("",@alp);
my $alp_regular = join("",@alp_out);
##right
chomp $arp;
my $arp_length = length $arp;
my @arp = split "",$arp;
my @arp_out;
for (my $i = 0; $i <= $#arp; $i ++){
    if (exists $dege{$arp[$i]}){
        push @arp_out,$dege_out{$arp[$i]};
        $arp[$i] = $dege{$arp[$i]};
    }else{
        push  @arp_out,$arp[$i];
    }
}
$arp = join("",@arp);
my $arp_regular = join("",@arp_out);

#get input handle
my $ass_r1_strings = gzip_support($r1_file,'input',3);
my $ass_r2_strings = gzip_support($r2_file,'input',3);
my $ass_r1_handle = new IO::File "$ass_r1_strings" or die "$!";
my $ass_r2_handle = new IO::File "$ass_r2_strings" or die "$!";

#get output handle
my $r1_left_out_strings = gzip_support('read-tag_left_R1.fastq.gz','output',4);
my $r2_left_out_strings = gzip_support('read-tag_left_R2.fastq.gz','output',4);
my $r1_left_out_handle = new IO::File "$r1_left_out_strings" or die "$!";
my $r2_left_out_handle = new IO::File "$r2_left_out_strings" or die "$!";

my $r1_right_out_strings = gzip_support('read-tag_right_R1.fastq.gz','output',4);
my $r2_right_out_strings = gzip_support('read-tag_right_R2.fastq.gz','output',4);
my $r1_right_out_handle = new IO::File "$r1_right_out_strings" or die "$!";
my $r2_right_out_handle = new IO::File "$r2_right_out_strings" or die "$!";

my $r1_fail_strings = gzip_support('unmatch_read-tag_R1.fastq.gz','output',2);
my $r2_fail_strings = gzip_support('unmatch_read-tag_R2.fastq.gz','output',2);
my $r1_fail_handle = new IO::File "$r1_fail_strings" or die "$!";
my $r2_fail_handle = new IO::File "$r2_fail_strings" or die "$!";
my $log_handle = new IO::File ">$log_file" or die "$!";

#defined stat variable
my $assemble_total_num = 0;
my $left_total_num = 0;
my $right_total_num = 0;
my %pos_zero = (
    'left' => 0,
    'right' => 0,
);
my %pos_one = (
    'left' => 0,
    'right' => 0,
);
my %pos_two = (
    'left' => 0,
    'right' => 0,
);
my $correct_total_num = 0;
my $correct_left_pass_num = 0;
my $correct_right_pass_num = 0;
my $correct_other_fail_num = 0;

#go
while (<$ass_r1_handle>){
    my $head = $_;
    my $seq = <$ass_r1_handle>;
    my $middle = <$ass_r1_handle>;
    my $qual = <$ass_r1_handle>;
    my $head2 = <$ass_r2_handle>;
    my $seq2 = <$ass_r2_handle>;
    my $middle2 = <$ass_r2_handle>;
    my $qual2 = <$ass_r2_handle>;
    chomp $seq2;
    chomp $qual2;
    my @SeqThread = SeqThread($seq2);
    if ($SeqThread[0] == 0){
        if (length $SeqThread[1] ne length $qual2){
            my $qual_number = -(length $SeqThread[1]);
            $qual2 = substr $qual2,$qual_number;
        }
        if ($SeqThread[2] eq 'left'){
            $left_total_num ++;
            print $r1_left_out_handle "$head$seq$middle$qual";
            print $r2_left_out_handle "$head$SeqThread[1]","\n+\n","$qual2\n";
        }else{
            $right_total_num ++;
            print $r1_right_out_handle "$head$seq$middle$qual";
            print $r2_right_out_handle "$head$SeqThread[1]","\n+\n","$qual2\n";
        }
    }else{
        print $r1_fail_handle "$head$seq$middle$qual";
        print $r2_fail_handle "$head2$seq2","\n+\n","$qual2\n";
    }
    $assemble_total_num ++;
}
$r1_left_out_handle -> close;
$r2_left_out_handle -> close;
$r1_right_out_handle -> close;
$r2_right_out_handle -> close;
$r1_fail_handle -> close;
$r2_fail_handle -> close;

#print stat
my $successful_seq = $left_total_num + $right_total_num;
my $ratio_successful = 100 * $successful_seq / $assemble_total_num;
printf $log_handle "Assemble Library Total Seq: %d, Split Successful Seq: %d,pr: %.2f%%.\n",$assemble_total_num,$successful_seq,$ratio_successful;
my $ratio_left = 100 * $left_total_num / $assemble_total_num;
my $ratio_right = 100 * $right_total_num / $assemble_total_num;
printf $log_handle "Assemble Library Split Successful Stat, LEFT: %d,pr: %.2f%%. RIGHT: %d,pr: %.2f%%.\n",$left_total_num,$ratio_left,$right_total_num,$ratio_right;
##Pos stat
my $ratio_left_zero = 100 * $pos_zero{'left'} / $left_total_num;
my $ratio_left_one = 100 * $pos_one{'left'} / $left_total_num;
my $ratio_left_two = 100 * $pos_two{'left'} / $left_total_num;
printf $log_handle "Assemble Library LEFT Pos Stat, Zero: %d,pr: %.2f%%. One: %d,pr: %.2f%%. Two: %d,pr: %.2f%%.\n",$pos_zero{'left'},$ratio_left_zero,$pos_one{'left'},$ratio_left_one,$pos_two{'left'},$ratio_left_two;
my $ratio_right_zero = 100 * $pos_zero{'right'} / $right_total_num;
my $ratio_right_one = 100 * $pos_one{'right'} / $right_total_num;
my $ratio_right_two = 100 * $pos_two{'right'} / $right_total_num;
printf $log_handle "Assemble Library RIGHT Pos Stat, Zero: %d,pr: %.2f%%. One: %d,pr: %.2f%%. Two: %d,pr: %.2f%%.\n",$pos_zero{'right'},$ratio_right_zero,$pos_one{'right'},$ratio_right_one,$pos_two{'right'},$ratio_right_two;
##Correct stat
my $ratio_other_fail = 100 * $correct_other_fail_num / $correct_total_num;
printf $log_handle "Assemble Library Correct Stat, Correct total: %d,Unable to determine: %d,pr: %.2f%%.\n",$correct_total_num,$correct_other_fail_num,$ratio_other_fail;
my $ratio_left_pass = 100 * $correct_left_pass_num / $correct_total_num;
my $ratio_right_pass = 100 * $correct_right_pass_num / $correct_total_num;
printf $log_handle "Assemble Library LEFT and RIGHT Correct Stat, LEFT Correct Pass: %d,pr: %.2f%% ,RIGHT Correct Pass: %d,pr: %.2f%%.\n",$correct_left_pass_num,$ratio_left_pass,$correct_right_pass_num,$ratio_right_pass;

$log_handle -> close;

#sub
sub SeqThread {
    my $in_seq = shift;
    #limit postion: 0-2bp
    #return 0 is success, return 1 is fail
    if ($in_seq =~ /^($alp)/ or $in_seq =~ /^($arp)/){
        if ($1 eq $alp){
            $pos_zero{'left'} ++;
            return 0,$in_seq,'left';
        }else{
            $pos_zero{'right'} ++;
            return 0,$in_seq,'right';
        }
    }elsif ($in_seq =~ /\w($alp)(.+)$/ or $in_seq =~ /\w($arp)(.+)$/){
        my $ture_seq = $1.$2;
        if ($1 eq $alp){
            $pos_one{'left'} ++;
            return 0,$ture_seq,'left';
        }else{
            $pos_one{'right'} ++;
            return 0,$ture_seq,'right';
        }
    }elsif ($in_seq =~ /\w\w($alp)(.+)$/ or $in_seq =~ /\w\w($arp)(.+)$/){
        my $ture_seq = $1.$2;
        if ($1 eq $alp){
            $pos_two{'left'} ++;
            return 0,$ture_seq,'left';
        }else{
            $pos_two{'right'} ++;
            return 0,$ture_seq,'right';
        }
    }else{
        $correct_total_num ++;
        my %distance_left = ();
        my %distance_right = ();
        for (my $i = 0;$i <= 2;$i ++){
            my $source_seq = substr $in_seq,$i,20;
            my $distance = HanMingDistanceLeft($source_seq);
            $distance_left{$i} = $distance if $distance <= $max_distance;
            $distance = HanMingDistanceRight($source_seq);
            $distance_right{$i} = $distance if $distance <= $max_distance;
        }
        #left correct
        if (scalar keys %distance_left > 0){
            my $min_values = min(values %distance_left);
            my @min_number = grep { $distance_left{$_} == $min_values } keys %distance_left;
            if (scalar @min_number == 1){
                my $ture_seq_begin = $min_number[0] + 20;
                my $ture_seq = substr $in_seq,$ture_seq_begin;
                $ture_seq = $alp_regular.$ture_seq;
                $correct_left_pass_num ++;
                return 0,$ture_seq,'left';
            }else{
                my $ture_seq_begin = 2 + 20;
                my $ture_seq = substr $in_seq,$ture_seq_begin;
                $ture_seq = $alp_regular.$ture_seq;
                $correct_left_pass_num ++;
                return 0,$ture_seq,'left';
            }
        }
        #right correct
        if (scalar keys %distance_right > 0){
            my $min_values = min(values %distance_right);
            my @min_number = grep { $distance_right{$_} == $min_values } keys %distance_right;
            if (scalar @min_number == 1){
                my $ture_seq_begin = $min_number[0] + 20;
                my $ture_seq = substr $in_seq,$ture_seq_begin;
                $ture_seq = $arp_regular.$ture_seq;
                $correct_right_pass_num ++;
                return 0,$ture_seq,'right';
            }else{
                my $ture_seq_begin = 2 + 20;
                my $ture_seq = substr $in_seq,$ture_seq_begin;
                $ture_seq = $arp_regular.$ture_seq;
                $correct_right_pass_num ++;
                return 0,$ture_seq,'right';
            }
        }
        $correct_other_fail_num ++;
        return 1,undef,undef;
    }
}

sub HanMingDistanceLeft {
    my $source_seq = shift;
    my $distance = 0;
    for (my $i = 0;$i <= $alp_length - 1; $i ++){
        $distance ++ if substr($source_seq,$i,1) !~ /${alp[$i]}/;
        return 100 if $distance > $max_distance;
    }
    return $distance;
}

sub HanMingDistanceRight {
    my $source_seq = shift;
    my $distance = 0;
    for (my $i = 0;$i <= $arp_length - 1; $i ++){
        $distance ++ if substr($source_seq,$i,1) !~ /${arp[$i]}/;
        return 100 if $distance > $max_distance;
    }
    return $distance;
}
