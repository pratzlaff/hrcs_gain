#! /usr/bin/perl -w
use strict;

die "this program is wrong!";

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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> April 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Fit::LM;
use PDL::Graphics::PGPLOT;
use File::Temp;
use Math::Trig qw( pi );
use Lab;

use Getopt::Long;
my %default_opts = (
		    mincnts => 50,
		    sigfilter => 2, # +/- prelim fitted sigma to use in refined gaussian
		    sigcut => 2.354/2, # half of fwhm
		    absres => 1,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'info!', 'mincnts=i', 'sigfilter=f', 'sigcut=f', 'absres!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV or die "usage: $0 [options] binfile\n";

my $bin = shift;

my $it = Lab::InBinFile->new($bin) or die;

my @cols = qw( ytap ysubtap xtap xsubtap y1 y2 x1 x2 sum );

if ($opts{info}) {

  print join("\t", @cols),"\n";
  print join("\t", ('N')x@cols),"\n";

  while (my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) =
	 $it->next_subtap) {
    print join("\t",
	       $ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $y->sum,
	       ),"\n";
  }
}


# do a single fit of the first dataset
else {

  push @cols, qw( gnorm gsigma gmean resmean resrms );
  print join("\t", @cols),"\n";
  print join("\t", ('N')x@cols),"\n";

  while (
	 my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) = $it->next_subtap
	) {
    my $counts = $y->sum;

    next unless $counts >= $opts{mincnts};

    $_ = $_->float for $x, $y;

    my $mean = sum($x*$y)/$y->sum;

    my $params = zeroes(3);
    my $params_init = double($y->sum, 15, $mean);

    # get initial fit of single gaussian
    my $fit;
    ($fit, $params) = lmfit($x, $y, 1, \&one_normal, $params_init);
    my ($gnorm, $gsigma, $gmean) = map { $params->slice("($_)") } (0..2);

    # now use only pha channels within a couple sigma of the peak to
    # refine the fit for the peak
    my $index = which(abs($x-$gmean) <= $opts{sigfilter}*$gsigma);

    # fit the main peak again
    ($fit, $params) = lmfit($x->index($index), $y->index($index), 1, \&one_normal, $params_init);
    ($gnorm, $gsigma, $gmean) = $params->list;

    # recalculate fitted values for all channels
    $fit = $gnorm / $gsigma / sqrt(2*pi) * exp(-.5*(($x-$gmean)/$gsigma)**2);

    # residuals with cutoff to left of main peak
    my $res = $y-$fit;
    if ($opts{absres}) {
      $res = abs($res);
    } else {
      $res->where($res < 0) .= 0;
    }

    $index = which( ( $x > 0 ) & ( $x <= $gmean-$opts{sigcut}*$gsigma ) );

    $_ = $_->index($index) for $x, $y, $res;

    my $resmean = $res->sum ? sum($x * $res) / $res->sum : 0;
    my $resrms = $res->sum ? sqrt( sum( ($x - $resmean)**2 * $res) / $res->sum) : 0;

    print join("\t",
	       $ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $counts,
	       sprintf("%.2f",$gnorm),
	       sprintf("%.2f",$gsigma),
	       sprintf("%.2f",$gmean),
	       sprintf("%.2f",$resmean),
	       sprintf("%.2f",$resrms),
	      ),"\n";
  }
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
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
