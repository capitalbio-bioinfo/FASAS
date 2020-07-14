#!/usr/bin/env perl
=head
    FileName: Contig_Quality.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Contig_Quality.pl
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
    Input: Contigs and 16s rRNA sequence
    Output: Contigs with quality value
=cut

=what
    计划建立contig的质量值，序列比对contig后寻找每一个碱基比对的位置，然后按最高的质量值建立哈希
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $options_number = scalar @ARGV;
Getopt(
    'contig_file=s' => \$contig_file,
    'forward_seq=s' => \$forward_seq,
    'reverse_seq=s' => \$reverse_seq,
    #filter_contig
    'contig_length=i' => \$contig_lenhth,
    'n_filter'
    'help' => \$help,
) or die "$!";

if ($help or $options_number < 1){
    print "
    
";
    exit;
}

=step
    step1: read contig file and filter contig
    step2: read fastq and filter fastq
    step3: create temp_folder, print contig and fastq
    step4: build bowtie2 index, run bowtie2
    step5: read bowtie2 result, loop every contig build quality
    step6: report build stat, print log file
=cut

