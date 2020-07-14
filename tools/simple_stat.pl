#!/usr/bin/perl
use strict;
use warnings;
use IO::File;
use Getopt::Long;

my ($tn5_stat, $cycle_stat, $fasas_stat, $output_file, $help);

my $options_number = scalar @ARGV;
GetOptions(
    'tn5_stat=s' => \$tn5_stat,
    'cycle_stat=s' => \$cycle_stat,
    'fasas_stat=s' => \$fasas_stat,
    'output_file=s' => \$output_file,
    'help' => \$help,
) or die "$!";

if ($help or $options_number < 1){
    print "
    
";
    exit;
}

my %tn5 = %{read_table($tn5_stat)};
my %cycle = %{read_table($cycle_stat)};
my %fasas = %{read_table($fasas_stat)};

my @sample_name = sort {$a cmp $b} keys %fasas;

#format output
my @output_head = qw(Sample Direction RawReads CleanReads UMIReads UMINumber BinReads BinUMI ContigN50 GoodContigNumber);

my @all_result;
#LEFT
foreach my $every_sample (@sample_name){
    my $item = print_stat(\%tn5, \%cycle, \%fasas, $every_sample, 'left');
    push @all_result, $item;
}
#RIGHT
foreach my $every_sample (@sample_name){
    my $item = print_stat(\%tn5, \%cycle, \%fasas, $every_sample, 'right');
    push @all_result, $item;
}

my $output_handle = new IO::File ">$output_file" or die "$!";
print $output_handle join("\t",@output_head),"\n";
foreach my $every_result (@all_result){
    my @item = @{$every_result};
    print $output_handle join("\t",@item),"\n";
}
close $output_handle;

sub print_stat {
    my %tn5_ref = %{shift @_};
    my %cycle_ref = %{shift @_};
    my %fasas_ref = %{shift @_};
    my $sample_name = shift;
    my $direction = shift;

    my @result = ($sample_name);
    if ($direction =~ /left/){
        push @result, 'LEFT';
        push @result, $tn5_ref{$sample_name}{'Left RawReads'};
        push @result, $tn5_ref{$sample_name}{'Left CleanReads'};
        my $umi_reads = $fasas_ref{$sample_name}{'LEFTGoodSeqN.'} + $fasas_ref{$sample_name}{'LEFTBadSeqN.'};
        push @result, $umi_reads;
        push @result, $fasas_ref{$sample_name}{'LEFTUMIN.'};
        push @result, $tn5_ref{$sample_name}{'Last Reads'};
        push @result, $fasas_ref{$sample_name}{'16SNumber'};
        push @result, $fasas_ref{$sample_name}{'ContigN50'};
        push @result, $fasas_ref{$sample_name}{'GoodContigNumber'}
    }elsif ($direction =~ /right/){
        push @result, 'RIGHT';
        push @result, $tn5_ref{$sample_name}{'Right RawReads'};
        push @result, $tn5_ref{$sample_name}{'Right CleanReads'};
        my $umi_reads = $fasas_ref{$sample_name}{'RIGHTGoodSeqN.'} + $fasas_ref{$sample_name}{'RIGHTBadSeqN.'};
        push @result, $umi_reads;
        push @result, $fasas_ref{$sample_name}{'RIGHTUMIN.'};
        push @result, $tn5_ref{$sample_name}{'Last Reads'};
        push @result, $fasas_ref{$sample_name}{'16SNumber'};
        push @result, $fasas_ref{$sample_name}{'ContigN50'};
        push @result, $fasas_ref{$sample_name}{'GoodContigNumber'}
    }else{
        die "$!"
    }

    return \@result;
}

sub read_table {
    my $input = shift;
    my $input_handle = new IO::File "$input" or die "$!";
    my $head = <$input_handle>;
    chomp $head;
    my @head = split "\t",$head;
    for (my $i = 0; $i <= $#head; $i ++){
        if ($head[$i] =~ /length\s+>=/){
            $head[$i] = 'GoodContigNumber';
        }
    }
    my %data;
    while (<$input_handle>){
        chomp;
        my @line = split "\t",$_;
        for (my $i = 1; $i <= $#head; $i ++){
            $data{$line[0]}{$head[$i]} = $line[$i];
        }
    }
    close $input_handle;
    return \%data;
}
