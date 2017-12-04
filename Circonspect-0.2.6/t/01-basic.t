#! /usr/bin/env perl

use strict;
use Test::More tests => 2;

BEGIN {
  use_ok( 'Circonspect' );
}

my $fastafiles  = ['./t/data/100seq.fa.gz', './t/data/100seq.fa.gz'];
my $qualfiles   = undef;
my $discard     = 100;
my $trim        = 100;
my $outdir      = '.';
my $outprefix   = 'test';
my $outbasic    = undef;
my $outdetails  = '-e';
my $outstats    = '-t';
my $outace      = undef;
my $samplesize  = 50;
my $metapercent = undef;
my $samplecompo = undef;
my $repnumber   = 2;
my $mincoverage = undef;
my $seed        = undef;
my $asmprog     = 'Minimo';
my $minidentity = 98;
my $minoverlap  = 35;
my $mixed       = '-g';
my $cross       = '-x';
my $cplxthres   = 0;
my $cplxwindow  = undef;
ok( Circonspect::Circonspect( $fastafiles, $qualfiles, $discard, $trim,
    $outdir, $outprefix, $outbasic, $outdetails, $outstats, $outace, $samplesize,
    $metapercent, $samplecompo, $repnumber, $mincoverage, $seed, $asmprog,
    $minidentity, $minoverlap, $mixed, $cross, $cplxthres, $cplxwindow), 
    'basic run' );
unlink 'test-details.txt';
