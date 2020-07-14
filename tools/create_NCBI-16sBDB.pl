#!/usr/bin/perl
=head
    FileName: Create_NCBI-16sBDB.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Create_NCBI-16sBDB.pl
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
    Input: NCBI Blast 16S fasta, NCBI Taxonomy databsae and NCBI accession2taxid file
    Output: Clear NCBI Blast 16S database file
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;
use Data::Dumper;
use File::Path qw(make_path);
use List::Util qw(uniq min);
use apackages;

my ($fasta_ref,$accession2taxid_file,$taxdump_names,$taxdump_nodes,$taxdump_merged,$output_dir,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "fasta_ref=s" => \$fasta_ref,
    "accession2taxid=s" => \$accession2taxid_file,
    "taxdump_names=s" => \$taxdump_names,
    "taxdump_nodes=s" => \$taxdump_nodes,
    "taxdump_merged=s" => \$taxdump_merged,
    "output_dir=s" => \$output_dir,
) or die "$!";

#main
##check options;
if ($help or $options_number < 1){
    ###help
    print "NCBI Blast 16S Database Generator
Usage:
    perl $0 -f [input fasta] -a [ID conversion file] -o [output folder] --taxdump_names [] --taxdump_nodes [] --taxdump_merged []
Options:
    fasta_ref         strings    Sequence files, extracted from the BLAST database using the blastdbcmd command. [require]
    accession2taxid   strings    File containing NCBI login number starting with N and Taxonomy database ID. [require]
    taxdump_names     strings    Names.dmp of the NCBI Taxonomy database. [require]
    taxdump_nodes     strings    Nodes.dmp of the NCBI Taxonomy database. [require]
    taxdump_merged    strings    Merged.dmp of the NCBI Taxonomy database. [require]
    output_dir        strings    Output folder. [require]
    help              switch     Show this help message and exit. [off]
";
    exit;
}

echo_error("CreateNBD","$fasta_ref: not found!") if ! -e $fasta_ref;
echo_error("CreateNBD","$accession2taxid_file: not found!") if ! -e $accession2taxid_file;
echo_error("CreateNBD","$taxdump_names: not found!") if ! -e $taxdump_names;
echo_error("CreateNBD","$taxdump_nodes: not found!") if ! -e $taxdump_nodes;
echo_error("CreateNBD","$taxdump_merged: not found!") if ! -e $taxdump_merged;
make_path($output_dir) if ! -d $output_dir;

##run
my $hash_accession = build_ref($fasta_ref);
my $hash_step2 = search_taxid($accession2taxid_file, $hash_accession);
my $hash_step3 = find_taxonomy($taxdump_names,$taxdump_nodes,$taxdump_merged,$hash_step2);
my $hash_step4 = seven_level($hash_step3);
my ($hash_step5, $hash_comp) = delete_dup($hash_step4);

my $fasta_output = $output_dir.'/16sMicrobial.fa';
my $taxonomy_output = $output_dir.'/16sMicrobial_taxonomy.tsv';
my $comp_output = $output_dir.'/16sMicrobial_group.tsv';

my $fasta_handle = new IO::File ">$fasta_output" or echo_error("CreateNBD","$!");
my $taxonomy_handle = new IO::File ">$taxonomy_output" or echo_error("CreateNBD","$!");
my %hash_step5 = %{$hash_step5};
foreach my $every_accession (%hash_step5){
    print $fasta_handle ">$every_accession\n$hash_step5{$every_accession}{'seq'}\n";
    print $taxonomy_handle "$every_accession\n$hash_step5{$every_accession}{'taxonomy'}\n";
}
$fasta_handle -> close;
$taxonomy_handle -> close;

my $comp_handle = new IO::File ">$comp_output" or echo_error("CreateNBD","$!");
print $comp_handle "TaxonomyName\tRankLevel\tConstituent\n";
my %comp_table = %{$hash_comp};
foreach my $every_taxonomy (keys %comp_table){
    foreach my $every_level (keys %{$comp_table{$every_taxonomy}}){
        print $comp_handle "$every_taxonomy\t$every_level\t$comp_table{$every_taxonomy}{$every_level}\n";
    }
}
$comp_handle -> close;

#sub
sub build_ref {
    my $input_fasta = shift;
    my %accession;
    my $fast_format = guess_fast_format($input_fasta);
    if ($fast_format eq 'complex_fasta' or $fast_format eq 'simple_fasta'){
        my $input_handle = new IO::File "$input_fasta" or echo_error("Build_ref","$!");
        while (my ($id, $decs, $seq) = read_complex_fast($input_handle)){
            $decs = $id.$decs;
            my @decs = split />/,$decs;
            foreach my $every_decs (@decs){
                next if $every_decs !~ /\|/;
                my @strings = split /\|/,$every_decs;
                my $name = '';
                if ($strings[4] =~ /\s+?(.*)\s+strain/){
                    $name = $1;
                }else{
                    my @name = split /\s+/,$strings[4];
                    $name = join(' ',$name[0],$name[1]);
                }
                $accession{$strings[3]}{'name'} = $name;
                $accession{$strings[3]}{'seq'} = $seq;
            }
            last if eof($input_handle);
        }
        $input_handle -> close;
    }else{
        echo_error("Build_ref","$fast_format: Input file is not fasta format!");
    }
    return \%accession;
}

sub search_taxid {
    my $input_a2tid = shift;
    my $hash_a = shift;
    my %hash = %{ $hash_a };
    my %a2tid;
    my $input_handle = new IO::File "$input_a2tid" or echo_error("Search_taxid","$!");
    my $dbug_handle = new IO::File ">search_nulltaxid.log" or echo_error("Search_taxid","$!");
    while (<$input_handle>){
        chomp;
        my @line = split /\t/,$_;
        $a2tid{$line[1]} = $line[2];
    }
    $input_handle -> close;
    foreach my $every_accession (keys %hash){
        if (exists $a2tid{$every_accession}){
            $hash{$every_accession}{'taxid'} = $a2tid{$every_accession};
        }else{
            print $dbug_handle "$every_accession\n";
        }
    }
    $dbug_handle -> close;
    return \%hash;
}

sub find_taxonomy {
    my $names = shift;
    my $nodes = shift;
    my $merged = shift;
    my $hash = shift;
    my %hash = %{ $hash };
    my $names_handle = new IO::File "$names" or echo_error("Find_taxonomy","$!");
    my $nodes_handle = new IO::File "$nodes" or echo_error("Find_taxonomy","$!");
    my $merged_handle = new IO::File "$merged" or echo_error("Find_taxonomy","$!");
    my %names;
    while (<$names_handle>){
        chomp;
        my @line = split /\t\|\t/,$_;
        $line[3] =~ s/\t\|// if $line[3] =~ /\t\|/;
        next if ($line[3] ne 'scientific name');
        $names{$line[0]}{'name'} = $line[1];
    }
    $names_handle -> close;
    while (<$nodes_handle>){
        chomp;
        my @line = split /\t\|\t/,$_;
        $names{$line[0]}{'pid'} = $line[1];
        $names{$line[0]}{'rank'} = $line[2];
    }
    $nodes_handle -> close;
    my %merged;
    while (<$merged_handle>){
        chomp;
        my @line = split /\t\|\t/,$_;
        $line[1] =~ s/\t\|// if $line[1] =~ /\t\|/;
        $merged{$line[0]} = $line[1];
    }
    $merged_handle -> close;
    foreach my $every_accession (keys %hash){
        my @taxonomy = ();
        my @rank = ();
        my $taxid = $hash{$every_accession}{'taxid'};
        REIF:
        if (exists $merged{$taxid}){
            $taxid = $merged{$taxid};
            $hash{$every_accession}{'taxid'} = $merged{$taxid};
        }
        if (exists $names{$taxid}){
            unshift @taxonomy,$names{$taxid}{'name'};
            unshift @rank,$names{$taxid}{'rank'};
            $taxid = $names{$taxid}{'pid'};
            goto REIF;
        }
        $hash{$every_accession}{'taxonomy'} = join(';',@taxonomy);
        $hash{$every_accession}{'rank'} = join(';',@rank);
    }
    return \%hash;
}

sub seven_level {
    my $hash = shift;
    my %hash = %{ $hash };
    my %level = (
        'kingdom' => 1,
        'phylum' => 1,
        'class' => 1,
        'order' => 1,
        'family' => 1,
        'genus' => 1,
        'species' => 1,
    );
    my %second_level = (
        'superkingdom' => 1,
        'superphylum' => 1,
        'superclass' => 1,
        'superorder' => 1,
        'superfamily' => 1,
        'supergenus' => 1,
        'superspecies' => 1,
        'subkingdom' => 1,
        'subphylum' => 1,
        'subclass' => 1,
        'suborder' => 1,
        'subfamily' => 1,
        'subgenus' => 1,
        'subspecies' => 1,
    );
    my @level = qw(species genus family order class phylum kingdom);
    my @array = qw(superkingdom kingdom subkingdom superphylum phylum subphylum superclass class subclass superorder order suborder superfamily family subfamily supergenus genus subgenus superspecies species subspecies);
    foreach my $every_accession (keys %hash){
        my @taxonomy = split ";",$hash{$every_accession}{'taxonomy'};
        my @rank = split ";",$hash{$every_accession}{'rank'};
        my @delete_taxonomy = ();
        my @delete_rank = ();
        for (my $i = 0; $i <= $#rank; $i ++){
            if (exists $level{$rank[$i]}){
                push @delete_taxonomy,$taxonomy[$i];
                push @delete_rank,$rank[$i];
            }
        }
        next if scalar @delete_rank == 7;
        if (scalar @delete_rank > 7){
            echo_error("Seven_level","$every_accession\t$hash{$every_accession}{'rank'} : Rnak Number gt 7!");
        }else{
            my %result = ();
            for (my $i = 0; $i <= $#delete_rank; $i ++){
                $result{$delete_rank[$i]} = $delete_taxonomy[$i];
            }
            my %second_ref = ();
            for (my $i = 0; $i <= $#rank; $i ++){
                $second_ref{$rank[$i]} = $taxonomy[$i];
            }
            my @lack_level;
            foreach my $every_level (@level){
                my @sure = grep { $_ eq $every_level } @delete_rank;
                push @lack_level,$every_level if scalar @sure == 0;
            }
            foreach my $every_lack (@lack_level){
                my @second_rank = grep { $_ =~ /$every_lack/ } @rank;
                if (scalar @second_rank == 1){
                    $result{$second_rank[0]} = $second_ref{$second_rank[0]};
                }elsif (scalar @second_rank > 1){
                    echo_error("Seven_level","$every_accession\t$hash{$every_accession}{'rank'} : Two identical rank!");
                }else{
                    $result{$every_lack} = 'no taxonomy';
                }
            }
            my @last_taxonomy = ();
            my @last_rank = ();
            foreach my $every_array (@array){
                push @last_taxonomy,$result{$every_array} if exists $result{$every_array};
                push @last_rank,$every_array if exists $result{$every_array};
            }
            $hash{$every_accession}{'taxonomy'} = join(";",@last_taxonomy);
            $hash{$every_accession}{'rank'} = join(";",@last_rank);
        }
    }
    return \%hash;
}

sub delete_dup {
    my $hash = shift;
    my %hash = %{ $hash };
    open DBUG,">>delete_dup.log";
    #修改%hash的值，返回%hash
    #build seq hash
    my %accession;
    foreach my $every_accession (keys %hash){
        my $seq = $hash{$every_accession}{'seq'};
        if ($accession{$seq}){
            $accession{$seq} .= ';'.$every_accession;
        }else{
            $accession{$seq} = $every_accession;
        }
    }
    my %comparison_table = ();
    foreach my $every_ (keys %accession){
        my @all_accession = split ";",$accession{$every_};
        if (scalar @all_accession > 1){
            my @species = ();
            my %dup_taxonomy = ();
            foreach my $every_accession (@all_accession){
                my @line = split ";",$hash{$every_accession}{'taxonomy'};
                $dup_taxonomy{$every_accession} = $hash{$every_accession}{'taxonomy'};
                push @species,$line[-1];
            }
            #species
            my $uniq_species_number = () = uniq(@species);
            if ($uniq_species_number == 1){
                for (my $l = 1; $l <= $#all_accession; $l ++){
                    delete $hash{$all_accession[$l]};
                }
            }elsif ($uniq_species_number > 1){
                #build dup rank hash
                my %dup_rank = ();
                foreach my $every_dup_tax (keys %dup_taxonomy){
                    my @line = split ";",$dup_taxonomy{$every_dup_tax};
                    for (my $i = 0; $i<=$#line; $i ++){
                        if (exists $dup_rank{$i}{$line[$i]}){
                            $dup_rank{$i}{$line[$i]} ++;
                        }else{
                            $dup_rank{$i}{$line[$i]} = 1;
                        }
                    }
                }
                #find dup rank
                my @dup_number = ();
                foreach my $every_level (keys %dup_rank){
                    push @dup_number,$every_level if scalar keys %{$dup_rank{$every_level}} > 1;
                }
                #find the smallest rank level
                my $min_dup_number = min(@dup_number);
                #change frist taxonomy
                my $ok_accession = $all_accession[0];
                my $ok_taxonomy = $hash{$ok_accession}{'taxonomy'};
                my @ok_taxonomy = split ";",$ok_taxonomy;
                my @origin_taxonomy = @ok_taxonomy;
                for (my $s = $#ok_taxonomy; $s > 0; $s --){
                    if($s >= $min_dup_number){
                        my $d = $s - 1;
                        $ok_taxonomy[$s] = $ok_taxonomy[$d].' group';
                    }else{
                        last;
                    }
                }
                $hash{$ok_accession}{'taxonomy'} = join(";",@ok_taxonomy);
                #save comparison table
                print DBUG join(";",@ok_taxonomy),"\n";
                for (my $k = $min_dup_number; $k <= 6; $k ++){
                    my $c_name = $ok_taxonomy[$k];
                    my $s = $k + 1;
                    $comparison_table{$c_name}{$s} = $all_accession[0].'::'.$origin_taxonomy[$k];
                    for (my $i = 1; $i <= $#all_accession; $i ++){
                        my @now_taxonomy = split ";",$hash{$all_accession[$i]}{'taxonomy'};
                        $comparison_table{$c_name}{$s} .= ";".$all_accession[$i].'::'.$now_taxonomy[$k];
                    }
                }
                #delete other taxonomy
                for (my $i = 1; $i <= $#all_accession; $i ++){
                    delete $hash{$all_accession[$i]};
                }
            }
        }
    }
    close DBUG;
    return \%hash,\%comparison_table;
}

__END__
