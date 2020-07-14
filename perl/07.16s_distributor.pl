#!usr/bin/perl
=head
    FileName: 07.Determine_16s.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    07.Determine_16s.pl : Part of the FASAS basic analysis script
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
    Input: Linked-Tag library statistcs data and Read-Tag library
    Output: All 16s rRNA sequence
=cut

use strict;
use warnings;
use IO::File;
use Switch;
use POSIX qw(strftime);
use Getopt::Long;
use List::Util qw(min);
use apackages;

my ($left_input_umi2name,$right_input_umi2name,$left_input_umitable,$right_input_umitable,$left_input_r1fastq,$right_input_r1fastq,$left_input_r2fastq,$right_input_r2fastq,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "left_input_umi2name=s" => \$left_input_umi2name,
    "right_input_umi2name=s" => \$right_input_umi2name,
    "left_input_umitable=s" => \$left_input_umitable,
    "right_input_umitable=s" => \$right_input_umitable,
    "left_input_r1fastq=s" => \$left_input_r1fastq,
    "right_input_r1fastq=s" => \$right_input_r1fastq,
    "left_input_r2fastq=s" => \$left_input_r2fastq,
    "right_input_r2fastq=s" => \$right_input_r2fastq,
    "help!" => \$help,
) or echo_error("16sDistributor","Failed to get options!");

if ($help or $options_number < 7){
    print "16S Sequence Distribution Program
Usage:
Options:
    
";
    exit;
}

echo_info("16sDistributor","BEGIN");

#read step6 umi2seqname
#use 25 as the sequence number threshold
echo_info("16sDistributor","Read umi2name file");
my %LEFT_UMI2SeqName = %{ ReadTsv($left_input_umi2name,'seqname') };
my %RIGHT_UMI2SeqName = %{ ReadTsv($right_input_umi2name,'seqname') };

#read step4 UMI information
#UMI's sequence support number is greater than or equal to 3 and is in UMI2SeqName
echo_info("16sDistributor","Read umitable file");
my ($StatUMI_1to2,$StatUMI_1to2_Num) = ReadTsv($left_input_umitable,'statumi1');
#StatUMI_1to2 and StatUMI_2to1 is Two-dimensional hash
#StatUMI_1to2_Num and StatUMI_2to1_Num is umi 2 number
my %StatUMI_1to2 = %{$StatUMI_1to2};
my %StatUMI_1to2_Num = %{$StatUMI_1to2_Num};
my ($StatUMI_2to1,$StatUMI_2to1_Num) = ReadTsv($right_input_umitable,'statumi2');
my %StatUMI_2to1 = %{$StatUMI_2to1};
my %StatUMI_2to1_Num = %{$StatUMI_2to1_Num};

##Main
my $dbug_output = './dbug.log';
my $dbug_handle = new IO::File ">$dbug_output" or echo_error("16sDistributor","Failed to open DBUG handle");

my %left;
my $LEFT_ERROR_handle = new IO::File ">LEFT_ERROR.log" or echo_error("16sDistributor","Failed to open LEFT_ERROR handle");
foreach my $line_LEFT_UMI2SeqName (keys %LEFT_UMI2SeqName){
    my @search_result = ();

    print $dbug_handle "-----------------------------\n";
    print $dbug_handle "Tag: $line_LEFT_UMI2SeqName\n";

    $line_LEFT_UMI2SeqName = reverse_complement($line_LEFT_UMI2SeqName) if ! exists $StatUMI_1to2{$line_LEFT_UMI2SeqName};
    if (exists $StatUMI_1to2{$line_LEFT_UMI2SeqName}){

        print $dbug_handle "StatUMI_1to2 exists\n";

        my %StatUMI_1to2_line = %{$StatUMI_1to2{$line_LEFT_UMI2SeqName}};
        my $range = 3;
        if (scalar keys %StatUMI_1to2_line < $range){
            $range = scalar keys %StatUMI_1to2_line;
        }

        print $dbug_handle "range: $range goto for\n";
        print $dbug_handle "\tU1T\tU2T\tU1N\tU2N\tU2\n";

        my $counter = 0;
        foreach my $every_u2_key (sort {$StatUMI_1to2_line{$b} <=> $StatUMI_1to2_line{$a}} keys %StatUMI_1to2_line){

            print $dbug_handle "\t$StatUMI_1to2_Num{$line_LEFT_UMI2SeqName}\t";
            print $dbug_handle "$StatUMI_2to1_Num{$every_u2_key}\t";
            print $dbug_handle "$StatUMI_1to2_line{$every_u2_key}\t";
            print $dbug_handle "$StatUMI_2to1{$every_u2_key}{$line_LEFT_UMI2SeqName}\t";
            print $dbug_handle "$every_u2_key\n";

            if (CheckRange($StatUMI_1to2_Num{$line_LEFT_UMI2SeqName},$StatUMI_1to2_line{$every_u2_key}) == 0){
                if ($#search_result == -1){
                    print $dbug_handle "\tsearch_result -1: ok\n";
                    if (exists $StatUMI_2to1{$every_u2_key}){
                        my $Check2to1_result = Check2to1($StatUMI_2to1_Num{$every_u2_key},$StatUMI_1to2_Num{$line_LEFT_UMI2SeqName},$StatUMI_2to1{$every_u2_key}{$line_LEFT_UMI2SeqName},$StatUMI_1to2_line{$every_u2_key});
                        print $dbug_handle "\tCheck2to1_result: $Check2to1_result\n";
                        if ($Check2to1_result == 0){
                            push @search_result,$every_u2_key;
                        }elsif ($Check2to1_result == 1){
                            print $dbug_handle "LEFT: $line_LEFT_UMI2SeqName -- $every_u2_key 2to1 CHECK fail\n";
                        }
                    }
                }elsif ($#search_result == 0){
                    print $dbug_handle "\tsearch_result 0: ok\n";
                    #检查当前U2值是不是和第一个U2的值接近
                    my $min_value = int($StatUMI_1to2_line{$search_result[0]} - $StatUMI_1to2_line{$search_result[0]} * 0.2);
                    if ($StatUMI_1to2_line{$every_u2_key} >= $min_value){
                        if (exists $StatUMI_2to1{$every_u2_key}){
                            my $Check2to1_result = Check2to1($StatUMI_2to1_Num{$every_u2_key},$StatUMI_1to2_Num{$line_LEFT_UMI2SeqName},$StatUMI_2to1{$every_u2_key}{$line_LEFT_UMI2SeqName},$StatUMI_1to2_line{$every_u2_key});
                            print $dbug_handle "\tCheck2to1_result: $Check2to1_result\n";
                            if ($Check2to1_result == 0){
                                push @search_result,$every_u2_key;
                            }elsif ($Check2to1_result == 1){
                                print $dbug_handle "LEFT: $line_LEFT_UMI2SeqName -- $every_u2_key 2to1 CHECK fail\n";
                            }
                        }
                    }
                }elsif ($#search_result == 1){
                    print $dbug_handle "\tsearch_result 1: ok\n";
                    my $min_value = int($StatUMI_1to2_line{$search_result[1]} - $StatUMI_1to2_line{$search_result[1]} * 0.2);
                    if ($StatUMI_1to2_line{$every_u2_key} >= $min_value){
                        print $dbug_handle "\tsearch_result 1: end\n";
                        print $dbug_handle "LEFT: $line_LEFT_UMI2SeqName have three corresponding RIGHT UMI!\n";
                        next;
                    }
                }
            }else{
                print $dbug_handle "\tCheckRange: fail\n";
            }
            $counter ++;
            last if $counter >= $range;
        }
    }else{
        print $LEFT_ERROR_handle "LEFT: $line_LEFT_UMI2SeqName 1to2 ERROR\n";
        next;
    }
    if ($#search_result == -1){
        print $dbug_handle "Failed to find RIGHT and input is: $line_LEFT_UMI2SeqName\n";
        next;
    }
    print $dbug_handle "-----------------------------\n";
    $left{$line_LEFT_UMI2SeqName} = [@search_result];
}

$LEFT_ERROR_handle -> close;
$dbug_handle -> close;

#check dup map
my %dup_right;
foreach my $every_left (keys %left){
    #find dup right umi
    next if scalar @{$left{$every_left}} < 2;
    foreach my $every_right (@{$left{$every_left}}){
        #if every_right have two or more U1, delete this left, and warning;
        my $right_true_number = 0;
        if ($StatUMI_2to1_Num{$every_right}){
            foreach my $every_right_line (keys %{$StatUMI_2to1{$every_right}}){
                if (CheckRange($StatUMI_2to1_Num{$every_right},$StatUMI_2to1{$every_right}{$every_right_line}) == 0){
                    $right_true_number ++;
                }
            }
        }
        if ($right_true_number > 1){
            $dup_right{'left'}{$every_left} = 1;
            foreach my $every_left_value (@{$left{$every_left}}){
                $dup_right{'right'}{$every_left_value} = 1;
            }
        }
    }
}

#filter dup right
my $delete_dup_handle = new IO::File ">delete_dup.log" or echo_error("16sDistributor","$!");
foreach my $every_left (keys %left){
    if ($dup_right{'left'}{$every_left}){
        print $delete_dup_handle "delete $every_left\t";
        print $delete_dup_handle join ("\t",@{$left{$every_left}}),"\n";
        delete $left{$every_left};
        next;
    }
    foreach my $every_right (@{$left{$every_left}}){
        if ($dup_right{'right'}{$every_right}){
            print $delete_dup_handle "delete $every_left\t";
            print $delete_dup_handle join ("\t",@{$left{$every_left}}),"\n";
            delete $left{$every_left};
            last;
        }
    }
}
$delete_dup_handle -> close;

#format hash left
my %format_left;
my %format_right;
#hash format
#format_left: LEFT UMI => join(LU + All RU)
#format_right: RIGHT UMI => join(LU + All RU)
foreach my $every_left (keys %left){
    if (ref($left{$every_left}) eq 'ARRAY'){
        $format_left{$every_left} = join("_",$every_left,@{$left{$every_left}});
        foreach my $every_right (@{$left{$every_left}}){
            $format_right{$every_right} = join("_",$every_left,@{$left{$every_left}});
        }
    }
}

#print result
#print format: fastq format
#@[LEFT UMI]_[ARRAY RIGHT UMI]_[NUMBER]
our %number;
PrintFastq($left_input_r1fastq,$left_input_r2fastq,\%format_left,'l');
PrintFastq($right_input_r1fastq,$right_input_r2fastq,\%format_right,'r');

sub CheckRange {
    my $tag_number = shift;
    my $input_number = shift;
    switch($tag_number){
        case 1 {return 0}
        case [2..6] {if ($input_number > 1){return 0;}}
        case [7..10] {if ($input_number > 2){return 0;}}
        else {if ($input_number > int($tag_number * 0.3)){return 0;}}
    }
    return 1;
}

sub Check2to1 {
    my $StatUMI_2to1_Total = shift;
    my $StatUMI_1to2_Total = shift;
    my $StatUMI_2to1_Val = shift;
    my $StatUMI_1to2_Val = shift;
    #(12Val/12Total) * (21Val/21Total) > 0.09
    my $value = ($StatUMI_1to2_Val * $StatUMI_2to1_Val / $StatUMI_1to2_Total / $StatUMI_2to1_Total);
    $value >= 0.09 ? return 0 : return 1;
}

sub ReadTsv {
    my $file = shift;
    my $type = shift;
    if ($type eq 'statumi1'){
        my $statumi_handle = new IO::File "$file" or echo_error("16sDistributor","$!");
        my %statumi;
        my %statu1_num;
        #In the linked-tag library, we do not place any restrictions on LEFT UMI
        while (<$statumi_handle>){
            chomp;
            my @line = split "\t",$_;
            next if $line[1] < 3;
            $line[0] = reverse_complement($line[0]);
            next if ! exists $LEFT_UMI2SeqName{$line[0]};
            $statu1_num{$line[0]} = $line[1];
            for (my $i = 2;$i <= $#line;$i += 2){
                $statumi{$line[0]}{$line[$i]} = $line[$i+1] if exists $RIGHT_UMI2SeqName{$line[$i]};
            }
        }
        $statumi_handle -> close;
        return \%statumi,\%statu1_num;
    }elsif ($type eq 'statumi2'){
        my $statumi_handle = new IO::File "$file" or echo_error("16sDistributor","$!");
        my %statumi;
        my %statu2_num;
        while (<$statumi_handle>){
            chomp;
            my @line = split "\t",$_;
            next if ! exists $RIGHT_UMI2SeqName{$line[0]};
            $statu2_num{$line[0]} = $line[1];
            for (my $i = 2;$i <= $#line;$i += 2){
                $line[$i] = reverse_complement($line[$i]);
                $statumi{$line[0]}{$line[$i]} = $line[$i+1] if exists $LEFT_UMI2SeqName{$line[$i]};
            }
        }
        $statumi_handle -> close;
        return \%statumi,\%statu2_num;
    }elsif ($type eq 'seqname'){
        #hash structure
        #umi => join("\t",seqname)
        my $seqname_handle = new IO::File "$file" or echo_error("16sDistributor","$!");
        my %seqname;
        while (<$seqname_handle>){
            chomp;
            #number of "\t" => a array => scalar => $num
            my $num = () = $_ =~ /\t/g;
            #define an average coverage of 5 layers
            #define length of sequence is 150 bp
            #number of UMI sequence support >= 25
            next if $num < 25;
            $_ =~ /^(.+?)\t(.*)$/ ? $seqname{$1} = $2 : echo_error("16sDistributor","$file error: regex error!\n$_\n");
        }
        $seqname_handle -> close;
        return \%seqname;
    }
}

sub PrintFastq {
    my $pf_input_r1 = shift;
    my $pf_input_r2 = shift;
    my $pf_ref_hash = shift;
    my $pf_type = shift;
    my %pf_ref_hash = %{$pf_ref_hash};
    print scalar keys %pf_ref_hash;
    my $r1_strings = gzip_support($pf_input_r1,'input',2);
    my $r2_strings = gzip_support($pf_input_r2,'input',2);
    my $r1_handle = new IO::File "$r1_strings" or echo_error("16sDistributor","$!");
    my $r2_handle = new IO::File "$r2_strings" or echo_error("16sDistributor","$!");
    my $r1_output_handle;
    my $r2_output_handle;
    if ($pf_type eq 'l'){
        my $r1_output_strings = gzip_support('all_sequence_R1.fastq.gz','output',4);
        my $r2_output_strings = gzip_support('all_sequence_R2.fastq.gz','output',4);
        $r1_output_handle = new IO::File "$r1_output_strings" or echo_error("16sDistributor","$!");
        $r2_output_handle = new IO::File "$r2_output_strings" or echo_error("16sDistributor","$!");
    }elsif ($pf_type eq 'r'){
        my $r1_output_strings = gzip_support('all_sequence_R1.fastq.gz','append',4);
        my $r2_output_strings = gzip_support('all_sequence_R2.fastq.gz','append',4);
        $r1_output_handle = new IO::File "$r1_output_strings" or echo_error("16sDistributor","$!");
        $r2_output_handle = new IO::File "$r2_output_strings" or echo_error("16sDistributor","$!");
    }else{
        echo_error("16sDistributor","Unknown type: $pf_type");
    }

    my $bad_output = new IO::File ">distributor_bad.log" or echo_error("16sDistributor","$!");
    while (<$r1_handle>){
        chomp;
        my $header2 = <$r2_handle>;
        my $seq = <$r1_handle>;
        my $seq2 = <$r2_handle>;
        my $middle = <$r1_handle>;
        my $middle2 = <$r2_handle>;
        my $qual = <$r1_handle>;
        my $qual2 = <$r2_handle>;
        if ($header2 =~ /.+ ([ATCGN]+)/){
            my $pf_key_value = $1;
            if (exists $pf_ref_hash{$pf_key_value}){
                $number{$pf_ref_hash{$pf_key_value}} += 1;
                print $r1_output_handle "\@$pf_ref_hash{$pf_key_value}_$number{$pf_ref_hash{$pf_key_value}}\n";
                print $r1_output_handle "$seq$middle$qual";
                print $r2_output_handle "\@$pf_ref_hash{$pf_key_value}_$number{$pf_ref_hash{$pf_key_value}}\n";
                print $r2_output_handle "$seq2$middle2$qual2";
            }else{
                print $bad_output "$header2";
            }
        }else{
            echo_error("16sDistributor","$header2 no UMI info!");
        }
    }
    $r1_handle -> close;
    $r2_handle -> close;
    $r1_output_handle -> close;
    $r2_output_handle -> close;
    $bad_output -> close;
}
