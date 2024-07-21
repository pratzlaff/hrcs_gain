#! /usr/bin/perl -w
use strict;

use Lab;
use MyRDB;
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Fit::Polynomial;

my $b = 'p197061012';
#$b = 'p197061112';
$b = 'p197061208';
my $r1 = "test/$b.rdb";
my $r2 = "$Lab::ANALDIR/$b.rdb";

my ($n1, $sm1, $sr1, $smed1, $pm1, $pr1, $pmed1, $siqr1, $piqr1, $stm1, $ptm1) = MyRDB::rdb_cols($r1, 'n', 'smean', 'srms', 'smed', 'pmean', 'prms', 'pmed', 'siqr', 'piqr', 'stmean_5', 'ptmean_5');
my ($n2, $sm2, $sr2, $smed2, $pm2, $pr2, $pmed2) = MyRDB::rdb_cols($r2, 'n', 'smean', 'srms', 'smed', 'pmean', 'prms', 'pmed');

$_ = pdl $_ for $sm1, $sr1, $smed1, $sm2, $sr2, $smed2, $pm1, $pr1, $pmed1, $pm2, $pr2, $pmed2, $siqr1, $piqr1, $stm1, $ptm1, $n1, $n2;

printf "samp mean diff max: %f\n", abs($sm1-$sm2)->max;
printf "samp rms diff max: %f\n", abs($sr1-$sr2)->max;
printf "samp median diff max: %f\n", abs($smed1-$smed2)->max;

printf "pha mean diff max: %f\n", abs($pm1-$pm2)->max;
printf "pha rms diff max: %f\n", abs($pr1-$pr2)->max;
printf "pha median diff max: %f\n", abs($pmed1-$pmed2)->max;

dev '/xs', 3,3;
bin hist $sm1->where($sm1>0 & $n1 > 50);
bin hist $stm1->where($stm1>0 & $n1 > 50);
printf "samp trimmed mean std dev = %f\n", ($stm1->where($stm1>0&$n1>50)->stats)[1];
printf "samp mean std dev = %f\n", ($sm1->where($sm1>0&$n1>50)->stats)[1];

bin hist $sr1->where($sr1>0);
bin hist $siqr1->where($siqr1>0);

#my $i = which($sr1 > 0 & $siqr1 > 0 & $n1 > 400);
my $i = which($n1 > 500);
my $r = $sr1 / $siqr1;
bin hist $r->index($i);
points $stm1->index($i), $r->index($i);


points $sr1->index($i), $siqr1->index($i);
hold;
my ($yfit, $coeffs) = fitpoly1d($sr1->index($i), $siqr1->index($i),2);
line $sr1->index($i), $yfit;
release;
print $coeffs,"\n";
