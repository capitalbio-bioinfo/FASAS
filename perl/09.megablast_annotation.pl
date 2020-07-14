#!/usr/bin/perl
=head
    FileName: 09.Annotate_Tax.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    09.Annotate_Tax.pl : Part of the FASAS basic analysis script
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
    Input: Nogap contigs
    Output: Clear contigs and annotation result of megablast
=cut

use warnings;
use threads;
use Cwd qw(abs_path);
use Switch;
use IO::File;
use apackages;
use POSIX qw(strftime);
use File::Basename;
use File::Path qw(make_path remove_tree);
use Getopt::Long;

=description
    Run megablast annotations for one sample
=cut

=bug
    Chinese
    1. 在读取Blast比对结果时，由于Contig中可能的大片段Deletion导致这个Contig将会有两个主要比对 忽略
    2. 在这个脚本中应该输出给下游OTU分析管道的文件 解决
=cut

#defined options viailable and get default options.
my ($ref_fa,$ref_tax,$input_fasta,$Megablast_database,$cpu_number,$contig_length,$blast_result,$mini_identity,$blast_evalue,$help);
my $work_dir = './';
$mini_identity = 97;

#get options.
my $option_number = scalar @ARGV;
GetOptions (
    "input_fasta=s" => \$input_fasta,
    "work_dir=s" => \$work_dir,
    "ref_fa=s" => \$ref_fa,
    "ref_tax=s" => \$ref_tax,
    "blast_db=s" => \$Megablast_database,
    "blast_result" => \$blast_result,
    "cpu_number=i" => \$cpu_number,
    "blast_evalue=s" => \$blast_evalue,
    "contig_length=i" => \$contig_length,
    "mini_identity=f" => \$mini_identity,
    "help" => \$help,
) or echo_error("AnnotateTax","Failed to get options");

#help.
if ($help or $option_number == 0){
    print "Sequence Annotation by MegaBlast Program
Usage:
    perl $0 --work_dir [work_dir] --ref_fa [ref_fa] --ref_tax [ref_tax] --blast_db [ref_db] ...
Options:
    input_fasta     strings   input sequence file that fasta format [./contigs.fa]
    work_dir        path      Work folder [required]
    ref_fa          path      Reference FASTA file [required]
    ref_tax         path      Reference taxonomy file [required]
    blast_db        path      Megablast database [required]
    blast_evalue    number    Minimum BLAST E-value threshold [1e-20]
    blast_result    switch    Do not run the blast program, use existing results [off]
    mini_identity   number    Minimum similarity threshold, range 1 to 99 [97]
    help            switch    print this messages [off]
";
    exit;
}

#check options
echo_error("AnnotateTax","Input file $input_fasta cannot be read!") if ! -r $input_fasta;
echo_error("AnnotateTax","Ref_fa $ref_fa cannot be !") if ! -r $ref_fa;
echo_error("AnnotateTax","Ref_tax $ref_tax cannot be read!") if ! -r $ref_tax;

my $log_handle = new IO::File ">Ana_taxonomy.log" or die "$!";

#usearch11 for chimera inspection
my $script_abs_path = abs_path($0);
my $script_dir_path = dirname($script_abs_path);
my $bin_dir = $script_dir_path.'/../bin';
$ENV{'PATH'}.=":$bin_dir";
#check program
my $usearch11_check_stat = system("usearch11 --version >& /dev/null");
$usearch11_check_stat == 0 ? echo_info("AnnotateTax","Find usearch11!") : echo_error("AnnotateTax","Unexpected error, usearch11 not found");
my $blast_check_stat = system("blastn -version >& /dev/null");
$blast_check_stat == 0 ? echo_info("AnnotateTax","Find blastn!") : echo_error("AnnotateTax","Unexpected error, blastn not found");

#read ref tax
echo_info("AnnotateTax","==========READ REF TAX==========");
my %ref_tax = %{Read2Hash($ref_tax,"\t")};

#move to work dir
mkdir $work_dir if ! -d $work_dir;
chdir $work_dir;

#filter_contigs
echo_info("AnnotateTax","==========READ CONTIGS==========");

my $contigs_handle = new IO::File "$input_fasta" or echo_error("AnnotateTax","$!");
my (%contigs,%short_seq);
my $total_raw_contig_number = 0;
while (<$contigs_handle>){
    chomp;
    $_ =~ s/^>//;
    $_ = @{[$_ =~ /^(.+)\s/]}[0] if $_ =~ /\s/;
    my $contigs_seq = <$contigs_handle>;
    chomp $contigs_seq;
    $total_raw_contig_number ++;
    length $contigs_seq >= $contig_length ? ( $contigs{$_} = $contigs_seq ) : ( $short_seq{$_} = $contigs_seq );
}
$contigs_handle -> close;
echo_info("AnnotateTax","Total contig number is $total_raw_contig_number");
my $long_contigs_number = scalar keys %contigs;
echo_info("AnnotateTax","Contig of length gt $contig_length have $long_contigs_number");

echo_info("AnnotateTax","==========Running MegaBlast==========");

my $blast_result_address = $work_dir.'/contigs_blast6.fa';
echo_warn("AnnotateTax","Because the old alignment results are found, the program will skip the blast comparison process!") if -e $blast_result_address;
if (! $blast_result){
    #prepare a temp folder
    my $local_time = strftime "%Y%m%d_%H:%M:%S", localtime;
    my $tmpdir = '/tmp/AnaTax_'.$local_time.'_TmpDir';
    make_path($tmpdir) ? echo_info("AnnotateTax","Create TMP Dir $tmpdir") : echo_error("AnnotateTax","Create TMP Dir ERROR");

    my $available_mem_size = CheckMem();
    switch ($available_mem_size) {
        case -1 { echo_error("AnnotateTax","Failed to get available memory value, /proc/meminfo not exists") }
        case -2 { echo_error("AnnotateTax","Failed to get available memory value, unknown units") }
        case -3 { echo_error("AnnotateTax","Failed to get available memory value, no MemAvailable value") }
        case { $_[0] < 15 } { echo_warn("AnnotateTax","Currently available memory is small, running blast may take more time") }
    }
    my $available_cpu_number = CheckCPU() + 1;

    my @threads;
    #limit
    my $max_blast_pall = 10;
    my $max_blastn_threads = 14;
    #default value.
    my $threads_reads_number = 1000;
    my $every_mem_used = 3; #units is GB, when 60000 ref sequence

    #cpu check
    my $now_blastn_pall = int((scalar keys %contigs) / $threads_reads_number);
    $now_blastn_pall = 1 if $now_blastn_pall < 1;
    $now_blastn_pall = $available_cpu_number if $available_cpu_number < $now_blastn_pall;
    $now_blastn_pall = $max_blast_pall if $max_blast_pall < $now_blastn_pall;
    my $now_blastn_threads = int($available_cpu_number / $now_blastn_pall) + 1;
    $now_blastn_threads = $max_blastn_threads if $max_blastn_threads < $now_blastn_threads;

    #mem check
    my $total_mem = $now_blastn_pall * $every_mem_used;
    if ($available_cpu_number < $total_mem){
        while ($available_cpu_number < $total_mem){
            $now_blastn_pall --;
            $total_mem = $now_blastn_pall * $every_mem_used;
            echo_error("AnnotateTax","System memory is insufficient!") if $now_blastn_pall == 0;
        }
    }

    my $threads_remember = 1;
    my $threads_reads_remember = 1;
    my $thread_contig_name = $tmpdir.'/contig_'.$threads_remember.'.fa';
    my $thread_handle = new IO::File ">$thread_contig_name";
    foreach my $every_read (keys %contigs){
        print $thread_handle ">$every_read\n$contigs{$every_read}\n";
        if ($threads_reads_remember == $threads_reads_number){
            $threads_reads_remember = 0;
            $thread_handle -> close;
            push @threads,threads->create(\&RunMegaBlast,$thread_contig_name,$threads_remember,$now_blastn_threads);
            $threads_remember ++;
            if (scalar @threads == $now_blastn_pall){
                my $every_thread = shift @threads;
                $every_thread -> join();
            }
            #restart handle
            $thread_contig_name = $tmpdir.'/contig_'.$threads_remember.'.fa';
            $thread_handle = new IO::File ">$thread_contig_name" or echo_error("AnnotateTax","$!");
        }
        $threads_reads_remember ++;
    }
    if ( fileno($thread_handle) ){
        $thread_handle -> close;
        push @threads,threads->create(\&RunMegaBlast,$thread_contig_name,$threads_remember,$now_blastn_threads);
    }
    #join threads
    while (@threads){
        my $every_thread = shift @threads;
        $every_thread -> join();
    }

    #cat threads result
    echo_info("AnnotateTax","==========Binding MegaBlast Result==========");

    my $cat_stat = system("cat $tmpdir/*.blast6 > ./contigs_blast6.tsv");
    echo_error("AnnotateTax","Cat TMP Dir ERROR") if $cat_stat != 0;
    #remove tmp dir
    remove_tree($tmpdir) ? echo_info("AnnotateTax","Remove TMP Dir $tmpdir") : echo_error("AnnotateTax","Remove TMP Dir ERROR");
}

#read blast result
echo_info("AnnotateTax","==========Read MegaBlast Result==========");

my %blast_result = %{ReadBlastResult("contigs_blast6.tsv")};
my %top_identity_blast_result = %{TopScoreBR(\%blast_result)};
my %fasta_stand = %{ThreadStandInfo()};

echo_info("AnnotateTax","==========STAND FASTA==========");

my $stand_fasta_handle = new IO::File ">contigs_stand.fa" or die "$!";
foreach my $every_contigs_id (keys %contigs){
    $every_contigs_id =~ s/^>//;
    if (exists $fasta_stand{$every_contigs_id}){
        if (exists $top_identity_blast_result{$every_contigs_id}){
            my $ref_id = @{[keys %{$top_identity_blast_result{$every_contigs_id}}]}[0];
            my @top_result = @{$top_identity_blast_result{$every_contigs_id}{$ref_id}};
            my $begin_of_seq = 0;
            $top_result[4] > $top_result[5] ? ( $begin_of_seq = $top_result[5] ) : ( $begin_of_seq = $top_result[4] );
            my $sequence = substr $contigs{$every_contigs_id},$begin_of_seq,$top_result[1];
            if ($fasta_stand{$every_contigs_id} eq 'plus'){
                print $stand_fasta_handle ">$every_contigs_id\n$sequence\n";
            }elsif ($fasta_stand{$every_contigs_id} eq 'minus'){
                $sequence = reverse_complement($sequence);
                print $stand_fasta_handle ">$every_contigs_id\n$sequence\n";
            }else{
                print $log_handle "unkown stand: $every_contigs_id\n";
            }
        }else{
            print $log_handle "$every_contigs_id :no top identity blast result\n";
        }
    }else{
        print $log_handle "$every_contigs_id :has been deleted\n";
    }
}
$stand_fasta_handle -> close;

#usearch11 ref_chimer
echo_info("AnnotateTax","=========USEARCH11 UCHIME2_REF==========");

my $usearch_threads = 20;
$available_cpu_number = CheckCPU() + 1;
$usearch_threads = $available_cpu_number if $available_cpu_number < $usearch_threads;

my $usearch_command = 'usearch11 -uchime2_ref contigs_stand.fa -db '.$ref_fa.' -uchimeout contigs_stand.uchime -strand plus -mode high_confidence -threads '.$usearch_threads.' &> /dev/null';
my $usearch_status = system("$usearch_command");
$usearch_status == 0 ? echo_info("AnnotateTax","usearch run ok!") : echo_error("AnnotateTax","usearch return $usearch_status!");

#delete chimera
echo_info("AnnotateTax","=========DELETE CHIMERA==========");

my $usearch_result_handle = new IO::File "<contigs_stand.uchime" or die "$!";
my %usearch;
while (<$usearch_result_handle>){
    chomp;
    my @line = split "\t",$_;
    if ($line[2] eq '?' or $line[2] eq 'N'){
        $usearch{$line[0]} = 'N';
        next;
    }elsif ($line[2] eq 'Y'){
        $usearch{$line[0]} = 'Y';
        exists $top_identity_blast_result{$line[0]} ? $top_identity_blast_result{$line[0]}{'chimera'} = 'chimera' : print $log_handle "delete chimera: not found $line[0] in blast result!\n";
    }
}
$usearch_result_handle -> close;

echo_info("AnnotateTax","=========OUTPUT CONTIGS==========");

my $output_contig_handle = new IO::File ">contigs_clear.fa";
foreach my $every_contigs (keys %contigs){
    if (exists $usearch{$every_contigs}){
        print $output_contig_handle ">$every_contigs\n$contigs{$every_contigs}\n" if $usearch{$every_contigs} ne 'Y';
    }else{
        echo_warn("AnnotateTax","$every_contigs no usearch status!");
    }
}
$output_contig_handle -> close;

#anaylsis blast result.
echo_info("AnnotateTax","==========Annotating MegaBlast Result==========");

my (%fa_blast_status,%level_status);
#requires 90% of the sequences to be aligned
my $align_length_limit = $contig_length - $contig_length * 0.1;
foreach my $every_query_id (keys %top_identity_blast_result){
    next if exists $top_identity_blast_result{$every_query_id}{'chimera'};
    foreach my $every_ref_id (keys %{$blast_result{$every_query_id}}){
        my @one_result = @{$blast_result{$every_query_id}{$every_ref_id}};
        if ($one_result[1] > $align_length_limit){
                #cition:
                #Beye M , Fahsi N , Raoult D , et al. Careful use of 16S rRNA gene sequence similarity values for the identification of Mycobacterium species.[J]. New Microbes & New Infections, 2018, 22(C):24-29.]
                #Caporaso J G , Kuczynski J , Stombaugh J , et al. QIIME allows analysis of high-throughput community sequencing data[J]. Nature Methods, 2010.]
                #http://qiime.org/scripts/assign_taxonomy.html
                if ($one_result[8] < $blast_evalue and $one_result[0] >= $mini_identity){
                    $fa_blast_status{$every_query_id} = $every_ref_id;
                    $level_status{$every_query_id} = 'all';
                }else{
                    $fa_blast_status{$every_query_id} = $every_ref_id;
                    $level_status{$every_query_id} = 'unclassified';
                }
        }else{
            $fa_blast_status{$every_query_id} = 'tooshort';
        }
    }
}

#open output handle and output sequences that cannot be processed correctly and main comment file.
echo_info("AnnotateTax","==========PRINT Annotation Infomatation AND Unclassified Sequence==========");

my $short_contig_handle = new IO::File ">bad_alignment.fa" or echo_error("AnnotateTax","$!");
my $unclassified_contig_handle = new IO::File ">unclassified.fa" or echo_error("AnnotateTax","$!");
my $tax_output_handle = new IO::File ">primer_taxonomy.tsv" or echo_error("AnnotateTax","$!");

print $short_contig_handle "#This file records a sequence whose alignment length is less than the threshold.\n";
print $short_contig_handle "#possible reason: gap, higher mismatch or no mapping.\n";

print $unclassified_contig_handle "#This file records the sequence of 'unclassified'.\n";

foreach my $every_contigs_id (keys %contigs){
    if (exists $fa_blast_status{$every_contigs_id} and $fa_blast_status{$every_contigs_id} ne 'tooshort'){
        if (exists $level_status{$every_contigs_id} and $level_status{$every_contigs_id} eq 'all'){
            RefUndefined($fa_blast_status{$every_contigs_id}) if ! defined $ref_tax{$fa_blast_status{$every_contigs_id}};
            print $tax_output_handle "$every_contigs_id\t$fa_blast_status{$every_contigs_id}\t$ref_tax{$fa_blast_status{$every_contigs_id}}\n";
        }else{
            print $unclassified_contig_handle ">$every_contigs_id\n$contigs{$every_contigs_id}\n";
            print $tax_output_handle "$every_contigs_id\tN\tunclassified\n";
        }
    }else{
        if (! exists $top_identity_blast_result{$every_contigs_id}{'chimera'}){
            print $short_contig_handle ">$every_contigs_id\n$contigs{$every_contigs_id}\n";
        }
    }
}
$short_contig_handle -> close;
$tax_output_handle -> close;
$log_handle -> close;

sub ThreadStandInfo {
    my %seq_stand;
    #blast6 format
    foreach my $tsi_every_query_id (keys %top_identity_blast_result){
        foreach my $tsi_every_ref_id (keys %{$blast_result{$tsi_every_query_id}}){
            my @one_result = @{$blast_result{$tsi_every_query_id}{$tsi_every_ref_id}};
            my $query_stand = $one_result[4] - $one_result[5];
            my $ref_stand = $one_result[6] - $one_result[7];
            if ($query_stand < 0 and $ref_stand < 0){
                $seq_stand{$tsi_every_query_id} = 'plus';
            }elsif ($query_stand > 0 and $ref_stand > 0){
                $seq_stand{$tsi_every_query_id} = 'plus';
            }elsif ($query_stand == 0 or $ref_stand == 0){
                print $log_handle "thread stand: ? values in plus/minus $tsi_every_query_id\n";
            }else{
                $seq_stand{$tsi_every_query_id} = 'minus';
            }
        }
    }
    return \%seq_stand;
}

sub Read2Hash {
    #Read Taxonomy file
    my $input_file = shift;
    my $split_str = shift;
    my $input_handle = new IO::File "$input_file" or echo_error("AnnotateTax","$!");
    my %hash = ();
    while (<$input_handle>){
        chomp;
        next if $_ =~ /^#/;
        my @line = split "$split_str",$_;
        $hash{$line[0]} = $line[1];
    }
    $input_handle -> close;
    return \%hash;
}

sub ReadBlastResult {
    #Read Blast Result, deault file is 'contigs_blast6.tsv'
    my $file = shift;
    my $file_handle = new IO::File "$file" or echo_error("AnnotateTax","$file $!");
    my %hash = ();
    my %remember = ();
    while (<$file_handle>){
        chomp;
        my @line = split "\t",$_;
        my $query_id = shift @line;
        my $ref_id = shift @line;
        if (exists $remember{$query_id}{$ref_id}){
            if (@{$remember{$query_id}{$ref_id}}[0] > $line[1] and @{$remember{$query_id}{$ref_id}}[1] > $line[9]){
                print $log_handle "delete second blast result: $_\n";
                next;
            }
        }else{
            $remember{$query_id}{$ref_id} = [$line[1],$line[9]];
        }
        $hash{$query_id}{$ref_id} = [@line];
    }
    $file_handle -> close;
    #find unmatched sequence and output it to a file
    my $unmatch_handle = new IO::File ">./unmatch_sequence.fa" or echo_error("AnnotateTax",$!);
    print $unmatch_handle "#This file records the sequence that did not appear in the blast result.\n";
    foreach my $every_contig (keys %contigs){
        print $unmatch_handle "$every_contig\n$contigs{$every_contig}\n" and delete $contigs{$every_contig} if ! exists $hash{$every_contig};
    }
    $unmatch_handle -> close;
    return \%hash;
}

sub TopScoreBR {
    #Select the optimal alignment result based on Blast's Score and E-value
    my %hash = %{ shift @_ };
    foreach my $every_query_id (keys %hash){
        my $floor_evalue = 10;
        my $top_score = 0;
        my $good_ref_id = '';
        foreach my $every_ref_id (keys %{$hash{$every_query_id}}){
            my @one_result = @{$hash{$every_query_id}{$every_ref_id}};
            if ($one_result[9] > $top_score){
                $top_score = $one_result[9];
                $floor_evalue = $one_result[8];
                $good_ref_id = $every_ref_id;
            }elsif ($one_result[9] == $top_score and $one_result[8] < $floor_evalue){
                $top_score = $one_result[9];
                $floor_evalue = $one_result[8];
                $good_ref_id = $every_ref_id;
            }
        }
        foreach my $every_ref_id (keys %{$hash{$every_query_id}}){
            delete $hash{$every_query_id}{$every_ref_id} if $every_ref_id ne $good_ref_id;
        }
    }
    return \%hash;
}

sub RunMegaBlast {
    my $mega_input = shift;
    my $mega_number = shift;
    my $mega_threads = shift;
    my $mega_output = basename($mega_input,'.fa');
    $mega_output .= '.blast6';
    my $output_dir = dirname($mega_input);
    $mega_output = $output_dir.'/'.$mega_output;
    my $baslt_command = 'blastn -query '.$mega_input.' -db '.$Megablast_database.' -task megablast -use_index true -word_size 44 -num_alignments 5 -num_threads '.$mega_threads.' -outfmt 6 > '.$mega_output;
    my $blast_status = system("$baslt_command");
    echo_error("AnnotateTax","Thread$mega_number Megablast return $blast_status!") if $blast_status != 0;
    $blast_status == 0 ? return 0 : return 1;
}

sub RefUndefined {
    my $now_ref_id = shift;
    echo_error("AnnotateTax","Taxonomy Annotation file of Reference missing information
Ref Sequence ID: $now_ref_id
Please checking your Taxonomy Annotation file");
}
