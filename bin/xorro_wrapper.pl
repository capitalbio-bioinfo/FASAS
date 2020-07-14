#!/usr/bin/env perl -w

#----------------------------------------------------------------#

#XORRO -- XOR-based Read Overlapper

#This is the Perl-based wrapper for controlling the "xorro" 
#program. The XORRO software is made up of the C-based "xorro"
#application and this program "xorro_wrapper.pl".

#Please see README.txt for additional information.

#Software written by Russell J. Dickson

#----------------------------------------------------------------#


(my $help = <<HELP_DOC) =~ s/\t+//gm;

This is the perl wrapper for XORRO the XOR-based paired end Read Overlapper
Please see README.txt for more help

Arguments
"-i1" input fastq file name for file 1
"-i2" input fastq file name for file 2
"-a" minimum overlap length to be accepted (def. 10)
"-q" quality offset (def. 33)
"-r" maxmum iterations
"-n" value applied each round to starting read length in iterations (def. 27)
"-s" starting read length (def. reads file to find read length)
"-o" output filename (def. "overlap.fastq")
"-f" fully align reads with multiple iterations - F or T(def.)
"-k" keep the non-overlapped reads in files L_out.fastq and R_out.fastq - T or F(def.)
"-h" print the help information

To do a full iterative run with defaults:
xorro_wrapper.pl -i1 FILE1.fastq -i2 FILE2.fastq 

Output will be in files with the "overlap" prefix.
File from round 1: overlap1_XXX_10.fastq
File from round 2: overlap2_XXX_10.fastq
etc. where XXX is the number of nucleotides that were considered in each run.

HELP_DOC

use strict;
use Cwd 'abs_path';

#args
#mandatory
        #left file
	my $left_filename;
        #right file
	my $right_filename;

#optional
        #maximum number of iterations -r (default until finished
	my $max_iters = 999999;

	#value subtracted from the length in each successive round (advanced) -n (def = 64)
	my $add_const = 27;	

	#output filename root -o (default overlap)
	my $output_filename = "overlap.fastq";
	
	#overlap length -a (def = 10)
	my $accept_length = 10;

	#quality offset -q (def - 33)
	my $quality_offset = 33;

	#keep temp files -k (def = F)
	my $keep_temp = "F";
	
	#read_length for first input -s
	my $min_length = 70; #minimum overlap length (min 70)

	#maximum read length (end point)
	my $read_length;

	#TODO
	#-f T --- means do one quick runthrough (on by default?)
	my $full_run = "T";

if (@ARGV){
        my $i = 0;
        foreach my $item(@ARGV){
        $i++;
                if ($item eq "-i1"){     #input fastq file name
                        $left_filename = $ARGV[$i];
                }elsif ($item eq "-i2"){     #input fastq file name
                        $right_filename = $ARGV[$i];
                }elsif ($item eq "-a"){     #minimum overlap length to be accepted
                        $accept_length = $ARGV[$i];
                }elsif ($item eq "-q"){     #quality offset
                        $quality_offset = $ARGV[$i];
                }elsif ($item eq "-k"){     #keep temp files
                        $keep_temp = $ARGV[$i];
		}elsif ($item eq "-r"){     #maxmum iterations
                        $max_iters = $ARGV[$i];
		}elsif ($item eq "-n"){     #value applied each round in iterations
                        $add_const = $ARGV[$i];
		}elsif ($item eq "-s"){     #starting read length
                        $min_length = $ARGV[$i];
		}elsif ($item eq "-o"){     #output
                        $output_filename = $ARGV[$i];
		}elsif ($item eq "-f"){     #do a full run with multiple iterations
                        $full_run = $ARGV[$i];


                }elsif($item eq "-h"){  #print the help information
                        print STDERR "$help";
                        exit;
                }
        }
}else{
        print STDERR $help;
        exit;
}

#calculate read length if not defined
open FILE, $left_filename or die "Can't open file, $!\n";
my $line = <FILE>; #header
$line = <FILE>; #sequence
chomp $line;
my $actual_read_length = length($line);
if(!defined($read_length) || $actual_read_length < $read_length)
{
	$read_length = $actual_read_length;
}
#print "read_length = $read_length\n";
close FILE;

my $stop_length = $read_length;  #minimum overlap input length

open FILE, $right_filename or die "Can't open file, $!\n";
close FILE;


#the hard limit for shortest read length
if ($read_length < 70 || $stop_length < 70)
{
	die "Initial and minimum read length must be above 70\n";
}


my $main_loop_iter = 1;

#first run is on input files
#subsequent runs are on program-generated files
#print "$read_length\n";

#run once and exit
if($full_run ne "T")
{
	#print STDERR "xorro $left_filename $right_filename $min_length $accept_length $quality_offset > $output_filename\n"; #$main_loop_iter\_$read_length\_$accept_length\.fastq\n";
	#`./xorro $left_filename $right_filename $min_length $accept_length $quality_offset > $output_filename`; #$main_loop_iter\_$read_length\_$accept_length\.fastq`;

	exit;
}

my $install_location = abs_path($0);
$install_location =~ s/\/[^\/]*$//;

#print STDERR $install_location . "/xorro $left_filename $right_filename $min_length $accept_length $quality_offset > $output_filename\n"; #$main_loop_iter\_$min_length\_$accept_length\.fastq\n";
system $install_location . "/xorro $left_filename $right_filename $min_length $accept_length $quality_offset L_out.fastq R_out.fastq > $output_filename"; #$main_loop_iter\_$min_length\_$accept_length\.fastq";

$main_loop_iter++;
for(my $i = $min_length + $add_const; $i < $stop_length; $i += $add_const ) #start at second round in loop
{
	#move temporary output to input filenames
	#print "mv L_out.fastq L_in.fastq\n";
	#print "mv R_out.fastq R_in.fastq\n";
	`mv L_out.fastq L_in.fastq`;
	`mv R_out.fastq R_in.fastq`;

	#break out of main loop if the maximum number of rounds is reached
	if($main_loop_iter >= $max_iters)
	{
		#break out of loop
		last;
	}

	#run program - Updated to use absolute path
	my $install_location = abs_path($0);
	$install_location =~ s/\/[^\/]*$//;

	#`./xorro L_in.fastq R_in.fastq $i $accept_length $quality_offset L_out.fastq R_out.fastq >> $output_filename`; #$main_loop_iter\_$i\_$accept_length\.fastq`;	
    #print STDERR $install_location . "/xorro L_in.fastq R_in.fastq $i $accept_length $quality_offset L_out.fastq R_out.fastq >> $output_filename\n"; #$main_loop_iter\_$i\_$accept_length\.fastq\n";
	system $install_location . "/xorro L_in.fastq R_in.fastq $i $accept_length $quality_offset L_out.fastq R_out.fastq >> $output_filename";
	
	#remove old temporary files?
	#rename old temp files

	$main_loop_iter++;
}
if($main_loop_iter < $max_iters)
{
	`mv L_out.fastq L_in.fastq`;
	`mv R_out.fastq R_in.fastq`;
	
	my $install_location = abs_path($0);
	$install_location =~ s/\/[^\/]*$//;
	
    #print STDERR $install_location . "/xorro L_in.fastq R_in.fastq $stop_length $accept_length $quality_offset L_out.fastq R_out.fastq >> $output_filename\n"; #$main_loop_iter\_$stop_length\_$accept_length\.fastq\n";
	system $install_location . "/xorro L_in.fastq R_in.fastq $stop_length $accept_length $quality_offset L_out.fastq R_out.fastq >> $output_filename"; #$main_loop_iter\_$stop_length\_$accept_length\.fastq";
}

`rm L_in.fastq`;
`rm R_in.fastq`;

#remove temp files generated by xorro
if($keep_temp ne "T")
{
	`rm L_out.fastq`;
	`rm R_out.fastq`;
}
