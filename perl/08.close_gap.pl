#!/usr/bin/perl
=head
    FileName: 08.Fix_Contigs.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    08.Fix_Contigs.pl : Part of the FASAS basic analysis script
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
    Input: IDBA_ud assemble folder and 16s rRNA reference database
    Output: nogap contigs
=cut

=alignment
    bowtie2 end-to-end alignmnet
    all sequence is PLUS in the sam result file
=cut

=step_array
    only after IDBA_ud is completed
    dep IDBA_ud reuslt folder
=cut

=score_limit
    1. bowtie2 --score-min options
    2. bowtie2 default: L,-0.6,-0.6, score_min = int ( -0.6 - 0.6 * query_length )
    3. now default: L,-0.6,-0.4
=cut

=report_number
    1. bowtie2 -k options
    2. default value is 20
=cut

=bowtie2_sam
    YF:i: Read was filtered out?
    YI:Z: Summary of inputs to MAPQ calculation
    YM:i: Read was repetitive when aligned unpaired?
    YP:i: Read was repetitive when aligned paired?
    YT:Z: String representing alignment type
    YS:i: Score of other mate
    ZS:i: Pseudo-random seed
    XR:Z: Original read string
    XT:i: Time taken to align
    XD:i: DP problems
    XU:i: ungapped alignment
    YE:i: streak of failed DPs at end
    YL:i: longest streak of failed DPs
    YU:i: index of last succeeded DP
    YR:i: # redundant seed hits
    ZF:i: # FM Index ops
    ZI:i: # extend loop iters
    ZP:i: Score of best/second-best paired-end alignment
    ZU:i: Score of best/second-best unpaired alignment
    CP: Concordant; DP: Discordant; UP: Unpaired Mate; UU: Unpaired.
=cut

use strict;
use warnings;
use threads;
use apackages;
use POSIX qw(strftime);
use Fcntl qw(:flock);
use File::Path qw(make_path remove_tree);
use IO::File;
use Getopt::Long;
use List::Util qw(min max);
use Data::Dumper;

my ($idba_folder,$contigs,$db,$contig_limit,$max_threads,$output,$log,$help,$folder_log,$report_number,$score_min);

$idba_folder = './';
$folder_log = './log';
$contigs = './contigs.fa';
$db = '';
$contig_limit = 1200;
$max_threads = 20;
$report_number = 20;
$score_min = 'L,-0.6,-0.4';

#get options and help
my $option_number = scalar @ARGV;
GetOptions (
    "folder_log=s" => \$folder_log,
    "contigs=s" => \$contigs,
    "db=s" => \$db,
    "output=s" => \$output,
    "log=s" => \$log,
    "contig_limit=i" => \$contig_limit,
    "max_threads=i" => \$max_threads,
    "report_number=s" => \$report_number,
    "help!" => \$help,
) or echo_error("FixContigs","Failed to get options");

if ($help or $option_number == 0){
    print "
";
    exit;
}

#read idba_ud folder address
my $folder_log_handle = new IO::File "$folder_log" or echo_error("FixContigs","$!");
$idba_folder = <$folder_log_handle>;
chomp $idba_folder;

#read contigs
#%contigs: contig_name => contig_sequence
my $contigs_handle = new IO::File "$contigs" or echo_error("FixContigs","$!");

my $local_time = strftime "%Y%m%d_%H-%M-%S", localtime;
my $tmpdir = '/tmp/Second_Contig_'.$local_time.'_TmpDir';
make_path($tmpdir) ? echo_info("FixContigs","Create TMP Dir $tmpdir") : echo_error("FixContigs","Failed to Create TMP Dir");

#threads

my $trigger_counter = 0;

my %temp;
my @threads;
my $thread_number = 1;
while (<$contigs_handle>){
    chomp;
    my $contig_sequence = <$contigs_handle>;
    chomp $contig_sequence;
    echo_warn("FixContigs","Blank line:\n$_\n$contig_sequence") and next if $_ =~ /^$/;
    $temp{$_} = $contig_sequence;
    if (scalar keys %temp == 100){
        $thread_number = 1 if $thread_number > $max_threads;
        push @threads,threads -> create(\&bowtie2_contig,\%temp,$db,$tmpdir,$thread_number);
        $thread_number ++;
        %temp = ();
        if (scalar @threads == $max_threads){
            my $every_thread = shift @threads;
            $every_thread -> join();
        }
        $trigger_counter ++;
    }
}
$contigs_handle -> close;

if (scalar keys %temp >= 1){
    push @threads,threads -> create(\&bowtie2_contig,\%temp,$db,$tmpdir,$thread_number);
    $trigger_counter ++;
}

while (my $every_thread = shift @threads){
    $every_thread -> join();
}

echo_error("FixContigs","No Contig Input,exit!") if $trigger_counter == 0;

my $all_output = $tmpdir.'/*.fa';
my $all_warn = $tmpdir.'/*.txt';
my $cat_status = system("cat $all_output > $output");
$cat_status == 0 ? echo_info("FixContigs","Cat Result to $output") : echo_error("FixContigs","Failed to cat TMP Result!");
$cat_status = system("cat $all_warn > $log");
$cat_status == 0 ? echo_info("FixContigs","Cat LOG to $log") : echo_error("FixContigs","Failed to cat TMP LOG!");

#remove_tree($tmpdir) ? echo_info("FixContigs","Remove TMP folder $tmpdir") : echo_error("FixContigs","Remove TMP folder ERROR");
#remove_tree($idba_folder) ? echo_info("FixContigs","Remove IDBA folder $idba_folder") : echo_error("FixContigs","Remove IDBA folder ERROR");

#sub function
sub bowtie2_contig {
    my %data_hash = %{ shift @_ };
    my $bowtie2_db = shift;
    my $target_folder = shift;
    my $process_index = shift;

    my $output_file = $target_folder.'/output_file_'.$process_index.'.fa';
    my $warn_file = $target_folder.'/warn_file_'.$process_index.'.txt';
    my $output_handle = new IO::File ">>$output_file" or echo_error("FixContigs","$!");
    flock $output_handle,LOCK_EX;

    my %contigs_info;
    my %long_contigs;
    foreach my $every_contig (keys %data_hash){
        $long_contigs{$every_contig} = $data_hash{$every_contig} and next if length $data_hash{$every_contig} >= 1200;
        #only supports IDBA_ud results
        #sequence name format ex. CAGATCTGTTGAGACACGGC_42_1_1292_1292
        #umi_kmer_contig-number_contig-total-length
        if ($every_contig =~ />([a-zA-Z]+(?:_[a-zA-Z]+)*)_(\d+)_(\d+)_(\d+)_\d+/){
            $contigs_info{$every_contig}{'umi'} = $1;
            $contigs_info{$every_contig}{'kmer'} = $2;
            $contigs_info{$every_contig}{'contig_number'} = $3;
            $contigs_info{$every_contig}{'length'} = $4;
            $contigs_info{$every_contig}{'sequence'} = $data_hash{$every_contig};
        }else{
            echo_error("FixContigs","A regular expression error occurred while reading Contig!\n$every_contig");
        }
         my $contig_file = $idba_folder.'/'.$contigs_info{$every_contig}{'umi'}.'/contig-'.$contigs_info{$every_contig}{'kmer'}.'.fa';
        $contigs_info{$every_contig}{'file'} = $contig_file;
    }

    foreach my $name (keys %contigs_info){
        my %bowtie2_result;
        my %bowtie2_sequence;

        echo_error("FixContigs","Process: $process_index, $contigs_info{$name}{'file'} not found!") if ! -e $contigs_info{$name}{'file'};
        my $bowtie2_command = "bowtie2 -f $contigs_info{$name}{'file'} -x $bowtie2_db --mm --no-hd --no-unal --no-sq --quiet -k $report_number -p 3 --score-min $score_min";
        my $bowtie2_handle = new IO::File "$bowtie2_command |" or echo_error("FixContigs","Bowtie2 error, command:\n$bowtie2_command");

        while (<$bowtie2_handle>){
            chomp;
            my @line = split "\t",$_;
            $bowtie2_sequence{$line[0]} = $line[9];
            $bowtie2_result{$line[0]}{'direction'}{$line[2]} = $line[1];
            $bowtie2_result{$line[0]}{'postion'}{$line[2]} = $line[3];
            $bowtie2_result{$line[0]}{'cigar'}{$line[2]} = $line[5];
            $line[11] =~ /AS:i:(-?\d+)/ ? ( $bowtie2_result{$line[0]}{'score'}{$line[2]} = $1 ) : echo_error "Failed to get score!\nProcess: $process_index, File: $contigs_info{$name}{'file'}\n$_";
        }
        $bowtie2_handle -> close;

        #the best contig
        my $first_name = @{[sort keys %bowtie2_result]}[0];
        #idba_ud contig format
        #case1: no extend
        if (! defined $first_name or $first_name !~ /_0$/){
            print $output_handle "$name\n$contigs_info{$name}{'sequence'}\n";
            next;
        }
        #case2: extend
        #the best score
        my $max_first_score = max(values %{$bowtie2_result{$first_name}{'score'}});
        #the best ref
        my @max_score_ref = grep { $bowtie2_result{$first_name}{'score'}{$_} == $max_first_score } keys %{$bowtie2_result{$first_name}{'score'}};
        echo_error("FixContigs","Failed to get max_score_ref!\n$data_hash{$name}") if scalar @max_score_ref == 0;

        if (scalar @max_score_ref == 1){
            my $ref_result = filter_bowtie2_result(\%bowtie2_result,$max_score_ref[0],$name);
            my $extend_result = extend_contig($ref_result,\%bowtie2_sequence,$first_name,$warn_file);
            my $reduce_contig = reduce_nbase($extend_result,$contig_limit);
            #fix contig name information
            my $now_contig_length = length $reduce_contig;
            my $now_name = '';
            if ( $name =~ />([a-zA-Z]+(?:_[a-zA-Z]+)*)(_\d+_\d+_)\d+(_\d+)/ ){
                $now_name = '>'.$1.$2.$now_contig_length.$3;
            }
            $reduce_contig eq 'biggap' ? print $output_handle "$name\n$contigs_info{$name}{'sequence'}\n" : print $output_handle "$now_name\n$reduce_contig\n";
        }else{
            @max_score_ref = grep { $bowtie2_result{$first_name}{'direction'}{$_} == 0 or $bowtie2_result{$first_name}{'direction'}{$_} == 16 } @max_score_ref;
            my $ref_result = filter_bowtie2_result(\%bowtie2_result,$max_score_ref[0],$name);
            my $extend_result = extend_contig($ref_result,\%bowtie2_sequence,$first_name,$warn_file);
            my $reduce_contig = reduce_nbase($extend_result,$contig_limit);
            my $now_contig_length = length $reduce_contig;
            my $now_name = '';
            if ( $name =~ />([a-zA-Z]+(?:_[a-zA-Z]+)*)(_\d+_\d+_)\d+(_\d+)/ ){
                $now_name = '>'.$1.$2.$now_contig_length.$3;
            }
            $reduce_contig eq 'biggap' ? print $output_handle "$name\n$contigs_info{$name}{'sequence'}\n" : print $output_handle "$now_name\n$reduce_contig\n";
        }
    }

    foreach my $every_long (keys %long_contigs){
        print $output_handle "$every_long\n$long_contigs{$every_long}\n";
    }

    flock $output_handle,LOCK_UN;
    $output_handle -> close;
}

sub reduce_nbase {
    #$limit should not be too small as possible to avoid 'scalar @max_value > 1'
    my $result = shift;
    my $limit = shift;
    return 'biggap' if $result eq 'biggap';
    my @contig = split /N+/,$result;
    my @length = map { length $_ } @contig;
    if (max(@length) >= $limit){
        my @max_value = grep { length $_ == max(@length) } @contig;
        if (scalar @max_value == 0){
            my $max_value = max(@length);
            echo_error("FixContigs","Abnormal, Max length is $max_value but the value does not match in the \@contig\n$result");
        }elsif (scalar @max_value == 1){
            return $max_value[0];
        }else{
            my $ex = join('.+',@max_value);
            $result =~ /($ex)/;
            return $1;
        }
    }else{
        return $result;
    }
}

sub extend_contig {
    my %ref_hash = %{ shift @_ };
    my %sequence = %{ shift @_ };
    my $top_chioce = shift;
    my $log_file = shift;

    echo_error("FixContigs","Ref hash no info!") if scalar keys %ref_hash == 0;
    return $sequence{$top_chioce} if scalar keys %ref_hash == 1;

    my $warn_handle = new IO::File ">>$log_file" or echo_error("FixContigs","$!");
    flock $warn_handle,LOCK_EX;

    my %recontig = ();
    foreach my $every_contig (sort { length $sequence{$b} <=> length $sequence{$a} } keys %ref_hash){
        my @seq = split '',$sequence{$every_contig};
        my $postion = @{[values %{$ref_hash{$every_contig}{'postion'}}]}[0];
        my $cigar = @{[values %{$ref_hash{$every_contig}{'cigar'}}]}[0];
        my %cigar = %{ analysis_cigar($cigar,$postion) };
        my $seq_index = 0;
        foreach my $every_cigar (sort { $a <=> $b } keys %cigar){
            if (defined $recontig{$every_cigar}){
                if ($cigar{$every_cigar} eq 'M' or $cigar{$every_cigar} eq 'I'){
                    if ($recontig{$every_cigar} ne $seq[$seq_index]){
                        print $warn_handle "$top_chioce:\n$every_cigar\t$recontig{$every_cigar}\t$seq[$seq_index]\n";
                    }
                    $seq_index ++;
                }
            }else{
                if ($cigar{$every_cigar} eq 'M' or $cigar{$every_cigar} eq 'I'){
                    $recontig{$every_cigar} = $seq[$seq_index];
                    $seq_index ++;
                }elsif ($cigar{$every_cigar} eq 'D'){
                    $recontig{$every_cigar} = '-';
                }
            }
        }
    }
    my $recontig = format_contig(\%recontig);
    flock $warn_handle,LOCK_UN;
    $warn_handle -> close;

    $recontig =~ /N{20,}/i ? return 'biggap' : return $recontig;
}

sub filter_bowtie2_result {
    my $hash = shift;
    my %hash = %{ $hash };
    my $ref = shift;
    my $name = shift;
    my %filter_result;
    foreach my $every_contig (keys %hash){
        foreach my $every_ref (keys %{$hash{$every_contig}{'direction'}}){
            $filter_result{$every_contig}{'direction'}{$every_ref} = $hash{$every_contig}{'direction'}{$every_ref} if $every_ref eq $ref;
        }
        foreach my $every_ref (keys %{$hash{$every_contig}{'postion'}}){
            $filter_result{$every_contig}{'postion'}{$every_ref} = $hash{$every_contig}{'postion'}{$every_ref} if $every_ref eq $ref;
        }
        foreach my $every_ref (keys %{$hash{$every_contig}{'cigar'}}){
            $filter_result{$every_contig}{'cigar'}{$every_ref} = $hash{$every_contig}{'cigar'}{$every_ref} if $every_ref eq $ref;
        }
        foreach my $every_ref (keys %{$hash{$every_contig}{'score'}}){
            $filter_result{$every_contig}{'score'}{$every_ref} = $hash{$every_contig}{'score'}{$every_ref} if $every_ref eq $ref;
        }
    }
    return \%filter_result;
}

sub format_contig {
    my %hash = %{ shift @_ };
    my $contig;
    my $min_index = min(keys %hash);
    my $max_index = max(keys %hash);
    for (my $i = $min_index; $i <= $max_index; $i ++){
        defined $hash{$i} ? ( $contig .= $hash{$i} ) : ( $contig .= 'N' );
    }
    $contig =~ s/-//g;
    return $contig;
}

sub analysis_cigar {
    my $cigar = shift;
    my $postion = shift;
    my %cigar = ();
    if (my @cigar = $cigar =~ /(\d+)([A-Z])/g){
        for (my $i = 0; $i <= $#cigar; $i += 2){
            my $s = $i + 1;
            for (my $o = 0; $o <= $cigar[$i] - 1; $o ++){
                my $now_postion = $postion + $o;
                $cigar{$now_postion} = $cigar[$s];
            }
            $postion += $cigar[$i];
        }
    }
    return \%cigar;
}
