#! /usr/bin/perl -w
use strict;

=head1 NAME

template - A template for Perl programs.

=head1 SYNOPSIS

cp template newprog

=head1 DESCRIPTION

blah blah blah

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> May 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Fit::Polynomial;
use PDL::Graphics::PGPLOT;
use MyRDB;

use Getopt::Long;
my %default_opts = (
		    rdbdir => '/data/legs/rpete/cal/hrcs_gain/StatsData',
		    chip => 1,
		    type => 'samp',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'rdbdir=s', 'chip=i', 'type=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my @anodes = qw( B C O Al );
my @energies = (183, 277, 525, 1487);
my $energies = pdl \@energies;

my %factors = qw(
		 B  1.045
		 O  1.01
		 );

my %data;

for my $anode (@anodes) {

  my $mean = pdl [];
  my $err = $mean->copy;

  for my $chip (1..3) {

    my $file = "$opts{rdbdir}/$anode$chip.rdb";
    print STDERR "reading $file...";
    my ($m, $e) = rdb_data($file);
    print STDERR "done\n";

    $mean = append($mean, $m);
    $err = append($err, $e);
  }

  if (exists $factors{$anode}) {
    $_ *= $factors{$anode} for $err, $mean;
  }

  @{$data{$anode}}{'mean', 'err'} = ($mean, $err);
}

# divide all means by the average mean
my $meanmean = zeroes($data{B}{mean});
my $s2 = $meanmean->copy;
my $n = 0;
for (keys %data) {
  $s2 += $data{$_}{err} * $data{$_}{err};
  $n++;
  $meanmean += $data{$_}{mean};
}
$_ /= $n for $meanmean, $s2;
my $s = sqrt($s2);
$data{$_}{mean} /= $meanmean for keys %data;

my $means = cat( map $data{$_}{mean}, @anodes )->xchg(0,1)->sever;
my ($yfit, $coeffs) = fitpoly1d($energies->log10, $means, 2);

my $intercept = $coeffs->slice('(0)');
my $slope = $coeffs->slice('(1)');

dev 'linear4.ps/ps', 2,2;
bin hist($intercept), { xtitle => 'intercept' };
bin hist($slope), { xtitle => 'slope' };
points $meanmean, $intercept, { xtitle => 'mean mean', ytitle => 'intercept' };
points $meanmean, $slope, { xtitle => 'mean mean', ytitle => 'slope' };

for (keys %data) {
#  points
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

sub rdb_data {
  my $f = shift;

  my @cols = MyRDB::rdb_col_names($f) or die;

  my ($subs, $mean, $n, $rms, $err);

  my %cols = (
	      samp => { mean => 'stmean5', rms => 'srms' },
	      spi => { mean => 'spitmean5', rms => 'spirms' },
	      pha => { mean => 'ptmean5', rms => 'prms' },
	      );

  if (grep { $_ eq 'BG' } @cols) {
    my ($nnet, $bg);
    ($subs, $n, $nnet, $bg, $mean, $rms) = MyRDB::rdb_cols($f, qw( subs n nnet BG ), $cols{$opts{type}}{mean}, $cols{$opts{type}}{rms} );
    $_ = pdl $_ for $mean, $n, $rms, $nnet, $bg;
    $err = $rms * sqrt($n+$bg) / $nnet;
  }
  else {
    ($subs, $n, $mean, $rms) = MyRDB::rdb_cols($f, qw( subs n ), $cols{$opts{type}}{mean}, $cols{$opts{type}}{rms} );
    $_ = pdl $_ for $mean, $n, $rms;
    $err = $rms / sqrt($n);
  }

  my %err_factors = (
		     '3x1' => sqrt(3),
		     '3x3' => 3,
		      );

  while (my ($v, $f) = each %err_factors) {
    my $index = long [ grep { $subs->[$_] eq $v } 0..$#{$subs} ];
    (my $tmp = $err->index($index)) *= $f if $index->nelem;
  }

  return $mean, $err;

}
