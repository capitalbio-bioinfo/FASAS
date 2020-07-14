#!/usr/bin/perl
=head
    FilaName: Reverse_Complement.pl
    Auther: Ke Zhang
    Version: 1.0.0
    Date: 2019.10.10
=cut

=license
    Reverse_Complement.pl
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
    Input: Strings of Sequence
    Output: Reverse complement result
=cut

use strict;
use warnings;
use Getopt::Long;

my ($input_seq,$enter,$help);

my $options_number = scalar @ARGV;
GetOptions (
    "input_seq=s" => \$input_seq,
    "enter" => \$enter,
    "help" => \$help,
) or die "Failed to Get Options!\n";

if ($help or $options_number == 0){
    print "DNA Sequence Reverse Complement Program
Usage:
    perl $0 --input_seq|-i [input sequence]
Options:
    input_seq    strings    Input DNA sequence [required]
    enter        switch     Enter [off]
    help         switch     Print this message [off]

";
    exit;
}

chomp $input_seq;
die "RC: $input_seq doesn't seem to be DNA Sequence, because \"$1\" !\n" if $input_seq =~ /([^ATCGRYMKSWHBVDN]+)/i;
$input_seq = reverse $input_seq;
$input_seq =~ tr/ATCGRYMKSWHBVDN/TAGCYRKMWSDVBHN/;
$enter ? print $input_seq,"\n" : print $input_seq;
