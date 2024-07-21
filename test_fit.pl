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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> January 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Fit::LM qw( lmfit );
use PDL::Fit::Gaussian qw( fitgauss1d );
use PDL::IO::Misc;
use Math::Trig qw( pi );

use Getopt::Long;
my %default_opts = (
		    dev => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'dev=s', 'double!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my $txtfile = @ARGV ? shift : 'pha_histogram_example';

my ($x, $y) = rcols $txtfile;


my ($params_init, $fit, $params);
if ($opts{double}) {
  $params_init = float($y->sum/10, 45, 40, $y->sum, 15, 100);
  ($fit, $params) = lmfit($x, $y, 1, \&two_normals, $params_init);
}
else {
  $params_init = float($y->sum, 15, 100);
  ($fit, $params) = lmfit($x, $y, 1, \&one_normal, $params_init);
}

dev $opts{dev};
line $x, $y;
hold;
line $x, $fit, { linestyle => 2 };

print $params,"\n";

=begin comment

# compare with fitgauss1d fitted parameters
my ($cen, $pk, $fwhm, $back, $err, $gfit) = fitgauss1d($x, $y);
my $sigma = $fwhm / 2.3548;
my $norm = $pk * $sigma * sqrt(2*pi);

hold;
line $x, $gfit, { color => 'red' };
hold;

print "$norm  $sigma   $cen\n";

=cut

release;

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

sub two_normals {
  my ($x,$par,$ym,$dyda) = @_;

  my ($n1,$s1,$u1,$n2,$s2,$u2) = map { $par->slice("($_)") } (0..5);

  my $exp1 = 1 / $s1 / sqrt(2*pi) * exp(-0.5 * (($x-$u1)/$s1)**2 );
  my $exp2 = 1 / $s2 / sqrt(2*pi) * exp(-0.5 * (($x-$u2)/$s2)**2 );

  $ym .= $n1*$exp1 + $n2*$exp2;

  my (@dy) = map {$dyda->slice(",($_)")} (0..5);

  $dy[0] .= $exp1;                                       # partial wrt norm
  $dy[1] .= $n1 / $s1 * $exp1 * ((($x-$u1)/$s1)**2 - 1); # partial wrt sigma
  $dy[2] .= ($x-$u1) * $n1 / $s1 / $s1 * $exp1;          # partial wrt mean

  $dy[3] .= $exp2;
  $dy[4] .= $n2 / $s2 * $exp2 * ((($x-$u2)/$s2)**2 - 1);
  $dy[5] .= ($x-$u2) * $n2 / $s2 / $s2 * $exp2;
}

sub one_normal {
  my ($x,$par,$ym,$dyda) = @_;

  my ($n1,$s1,$u1) = map { $par->slice("($_)") } (0..2);

  my $exp1 = 1 / $s1 / sqrt(2*pi) * exp(-0.5 * (($x-$u1)/$s1)**2 );

  $ym .= $n1*$exp1;

  my (@dy) = map {$dyda->slice(",($_)")} (0..2);

  $dy[0] .= $exp1;                                       # partial wrt norm
  $dy[1] .= $n1 / $s1 * $exp1 * ((($x-$u1)/$s1)**2 - 1); # partial wrt sigma
  $dy[2] .= ($x-$u1) * $n1 / $s1 / $s1 * $exp1;          # partial wrt mean
}

