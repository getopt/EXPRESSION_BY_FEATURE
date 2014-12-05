#!/usr/bin/env perl
use strict;
use FileHandle;

my ($f1, $f2) = map { openfile($_); } @ARGV;
sub openfile {
   my $fn = shift;
   my $fh = new FileHandle "$fn", "r";
   die "no handle for $fn" unless $fh;
   return $fh;
}

my %selectids = map { chomp; ((split("\t", $_))[1], (split("\t", $_))[12]); } <$f1>; 
close($f1);

while(<$f2>){
    chomp;
    /gene_id "(.*)"; transcript_id/;
    die "unknown gene $1!\n" if !defined $selectids{$1};
    print $_, "gene_name \"", $selectids{$1}, "\" \n";
}
close $f2;

