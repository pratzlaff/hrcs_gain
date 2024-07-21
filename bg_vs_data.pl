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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> August 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Fit::Polynomial;
use Lab;
use Chandra::Tools::Common;

use Getopt::Long;
my %default_opts = (
		    dev => '/xs',
		    plotcoeffs => 1,
		    srccnts => 50,
		    bgcnts => 10,
		    rangeci => 0.75,
		    '3x3' => 1,
		    stat => 'mean',
		    type => 'samp',
		    bgsubtract => 1,
		    pdlfit => 1,
		    mask => 1,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'dev=s', 'bgfile=s', 'plotcoeffs!', 'swap!',
	   'type=s', '3x3!', 'srccnts=i', 'bgcnts=i', 'rangeci=f',
	   'stat=s', 'bgsubtract!', 'mask!', 'pdlfit!',
	   'one!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my @anodes = @ARGV ? @ARGV : @{+(Lab::lines())[0]};

$Lab::TYPE = $opts{type};

my $statfunc =
  $opts{stat} eq 'mean' ? \&Lab::hists_mean :
  $opts{stat} eq 'median' ? \&Lab::hists_median :
  die "statistic '$opts{stat}' is unrecognized";

dev $opts{dev}, $opts{one} ? () : (2, 4), { HardLW => 4 };

my $bg_hists = Lab::bg_hists($opts{bgfile}?($opts{bgfile}) : ());
$bg_hists = Lab::three_by_three($bg_hists) if $opts{'3x3'};

if ($opts{mask}) {
  (my $tmp = Lab::hists_indexND($bg_hists, scalar whichND(!Lab::active_mask_bg_noiffyv()))) .= 0;
}

#
# recreate the PHA values...BUT...rdl() pads arrays with zero to the
# size of the subtap with the most events, so anywhere else using this
# result must truncate the first dimension to values in $n
#
#my $pha = $hists->rld(sequence(long,256))->float;

#
# PHA histogram from all subtaps combined
#
# we'd like to use hist($pha, -.5, 255.5, 1), but once again rld()
# padding makes that incorrect (histogram has enormous erroneous spike
# at PHA=0)
#
# so instead,
#bin $hists->clump(1,2)->xchg(0,1)->sumover;
#
# also can be done with:
#   $hists->xchg(0,2)->sumover->sumover
#   $hists->xchg(0,2)->clump(2)->sumover
#

#
# histogram of events in each subtap
#
#bin hist $hists->sumover;

# image of events in each subtap
my $bg_n = $bg_hists->sumover;

# image of average in each subtap
my $bg_stats = $statfunc->($bg_hists);

my $histx_ul = $bg_hists->getdim(0) - .5;

# histogram of the background statistics
my ($bg_stats_histx, $bg_stats_histy)= hist $bg_stats->where($bg_n >= $opts{bgcnts}), -0.5, $histx_ul, 1;
#bin $bg_stats_histx, $bg_stats_histy;
my $bg_stats_to_use = stats_range($bg_stats_histx, $bg_stats_histy);

# average of median/means for every anode
my $stats_mean = zeroes($bg_stats);
my $stats_n = zeroes($bg_stats);

# hold fitted coefficients vs energy for each anode
my (@m, @b);

my @e = Lab::energies(@anodes);

my %fits;

for (0..$#anodes) {

  my $energy = $e[$_];
  my $anode = $anodes[$_];

  my (@x, @y, @n, @err);
  my $src_hists = Lab::src_hists($anode);
  $src_hists = Lab::three_by_three($src_hists) if $opts{'3x3'};

  if ($opts{mask}) {
    (my $tmp = Lab::hists_indexND($src_hists, scalar whichND(!Lab::active_mask_src_noiffyv()))) .= 0;
  }

  if ($opts{bgsubtract} and $anode eq 'B-Ka') {
    print STDERR "performing background subtraction for $anode\n";
    # FIXME : hard-coded exposure times
    my %exptimes = (
		    'B-Ka' => 18,
		    'O-Ka' => 3,
		    'C-Ka' => 2,
		    'Al-Ka' => 2,
		    );
    die $anode unless exists $exptimes{$anode};
    $src_hists -= $bg_hists * ($exptimes{$anode} / 66);
  }

  my $src_n = $src_hists->sumover;

  my $src_stats = $statfunc->($src_hists);

  # add the statistics from this source onto the mean statistics array
  my $i = whichND($src_n >= $opts{srccnts});
  my $tmp;
  ($tmp = $stats_mean->indexND($i)) += $src_stats->indexND($i);
  ($tmp = $stats_n->indexND($i))++;

  my $name = uc($opts{type});

  my $xtitle = "X-ray $opts{stat} $name";
  my $ytitle = "background $opts{stat} $name";

  $fits{$anode} = [ plot_stats(
			       $src_stats, $bg_stats,
			       $src_n, $bg_n,
			       $opts{srccnts}, $opts{bgcnts},
			       $xtitle, $ytitle, $anode,
			       1,
			       !$opts{one},
			      )
		  ];
  push @b, $fits{$anode}->[3][0];
  push @m, $fits{$anode}->[3][1];

}

my $i = whichND($stats_n > 0);
(my $tmp = $stats_mean->indexND($i)) /= $stats_n->indexND($i);
my $name = uc($opts{type});

my $xtitle = "average X-ray $opts{stat} $name";
my $ytitle = "background $opts{stat} $name";
$fits{Average} = [
		  plot_stats(
			     $stats_mean, $bg_stats,
			     $stats_n, $bg_n,
			     scalar(@anodes), $opts{bgcnts},
			     $xtitle, $ytitle, 'Average',
			     0,
			     !$opts{one},
			    )
		 ];

# Brad wanted everything on a single plot...
if ($opts{one}) {
  my @anodes = (@anodes, 'Average');
  my @x = map $fits{$_}[0], @anodes;
  my @y = map $fits{$_}[1], @anodes;
  my @yfit = map $fits{$_}[2], @anodes;
  my @coeffs = map $fits{$_}[3], @anodes;
  my @coeffs_errs = map $fits{$_}[4], @anodes;

  my $xmax = pdl(map $_->max, @x)->max;
  my $xmin = pdl(map $_->min, @x)->min;
  my $xdiff = $xmax - $xmin;
  $xmax += $xdiff * .05;
  $xmin -= $xdiff * .05;

  my $ymax = pdl(map $_->max, @y)->max;
  my $ymin = pdl(map $_->min, @y)->min;
  my $ydiff = $ymax - $ymin;
  $ymax += $ydiff * .05;
  $ymin -= $ydiff * .05;

  my @symbols = (2..@anodes, -8); # average gets filled circle
  my @colors = (2..@anodes, 1);   # average is black

  my $xtitle = "X-ray $opts{stat} $name";
  my $ytitle = "background $opts{stat} $name";
  swap($xtitle, $ytitle) if $opts{swap};

  for my $i (reverse(0..$#anodes)) { # plot combined first so axis titles are black
    points(
	   $x[$i],
	   $y[$i],
	   {
	    xrange => [$xmin, $xmax],
	    yrange => [$ymin, $ymax],
	    symbol => $symbols[$i],
	    color => $colors[$i],
	    xtitle => $xtitle,
	    ytitle => $ytitle
	   },
	  );
    hold;
    my $xline = zeroes($x[$i]->nelem+2);
    (my $tmp = $xline->slice("1:-2")) .= $x[$i];
    $xline->set(0, $xmin);
    $xline->set(-1, $xmax);
    my $yline = $coeffs[$i][0] + $coeffs[$i][1] * $xline;
    line($xline, $yline, { color => $colors[$i]});
    hold;
  }
#   legend(\@anodes
# 	 $xmin + $xdiff,
# 	 $ymax - $ydiff,
# 	 { color => \@colors, symbol => \@symbols },
# 	 );
  Chandra::Tools::Common::legend(
				 text => \@anodes,
				 points => \@symbols,
				 colors => \@colors,
				 lines => [(1)x@anodes],
				 left => 1,
				 );
  release;
}

if ($opts{plotcoeffs}) {

  my $e = (pdl(\@e) * 1e-3)->log10;

  points $e, pdl(\@m),
    {
     xtitle => 'energy',
     ytitle => 'slope',
     axis=>'logx',
     border => 1,
    };
  points $e, pdl(\@b),
    {
     xtitle => 'energy',
     ytitle => 'intercept',
     axis=>'logx',
     border => 1,
    };
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

sub stats_range {
  my ($histx, $histy) = @_;

  my ($stats_histx, $stats_histy) = ($histx, $histy);

  my $mean = Lab::hists_mean($histy, $histx);
  my $sdev = Lab::hists_prms($histy, $histx);

  return $histx->where(
		       ( $histx >= $mean-$sdev*sqrt(2)*erfi($opts{rangeci})) &
		       ( $histx <= $mean+$sdev*sqrt(2)*erfi($opts{rangeci}))
		       );
}

sub swap {
  my $tmp = $_[0];
  $_[0] = $_[1];
  $_[1] = $tmp;
}

sub plot_stats {

  my ($ind_stat, $dep_stat,
      $ind_n, $dep_n,
      $ind_cnts, $dep_cnts,
      $xtitle, $ytitle, $title,
      $weighted, $plotit,
     ) = @_;

  swap($ind_stat, $dep_stat);
  swap($ind_n, $dep_n);
  swap($ind_cnts, $dep_cnts);
  swap($xtitle, $ytitle);

  my $stats_to_use = stats_range(hist $ind_stat->where($ind_n >= $ind_cnts), -0.5, $histx_ul, 1);
  my (@x, @y, @n);
  for my $this_stat ($stats_to_use->list) {
    my $i = whichND(
		    ($ind_n>=$ind_cnts) &
		    ($dep_n>=$dep_cnts) &
		    (rint($ind_stat)==rint($this_stat)) &
		    1
		   );
    push @x, $this_stat;
    push @y, $dep_stat->indexND($i)->avg;
    push @n, $dep_n->indexND($i)->sum;
  }
  my ($x, $y, $n) = (pdl(\@x), pdl(\@y), pdl(\@n));

  if ($plotit) {
    points($x, $y,
	   {
	    xtitle => $xtitle,
	    ytitle => $ytitle,
	    title =>  $title,
	    border => 1,
	   }
	  );
    hold;
  }

  my ($yfit, @coeffs, @coeffs_err);
  if ($opts{pdlfit}) {
    my $w = $weighted ? $n : ones($n);
    my $coeffs;
    ($yfit, $coeffs) = fitpoly1d($x, $y, 2, { Weights => $w });
    @coeffs = $coeffs->list;
    @coeffs_err = (0, 0);
  }
  else {
    my ($b0, $b1, $s_b0, $s_b1);
    ($b0, $b1, $yfit, $s_b0, $s_b1) = linear_fit($x, $y);
    @coeffs = ($b0, $b1);
    @coeffs_err = ($s_b0, $s_b1);
  }

  if ($plotit) {
    line $x, $yfit;
    hold;
    text sprintf("\\fiy = %.3f (%.3f) x + %.3f (%.3f)",
		 $coeffs[1], $coeffs_err[1],
		 $coeffs[0], $coeffs_err[0]
		),
		  $x->max - ($x->max-$x->min)*0.1,
		    $y->min + ($y->max-$y->min)*0.1,
		      {
		       justification => 1 };
    release;
  }

  return $x, $y, $yfit, \@coeffs, \@coeffs_err;
}

# from stat 111 hw 11
sub linear_fit {
  my ($x, $y) = @_;
  my $xmean = $x->average;
  my $ymean = $y->average;

  my $n = $y->getdim(0);

  my $b1 = (($x-$xmean) * ($y-$ymean))->sumover / (($x-$xmean)**2)->sumover;
  my $b0 = $ymean - $b1 * $xmean;
  my $yfit = $b0 + $x * $b1;

  my $rss = (($yfit - $y)**2)->sumover;
  my $s = sqrt($rss / ($n-2));

  my $x2_sum = ($x * $x)->sumover;
  my $x_sum = $x->sumover;

  my $s_b0 = sqrt(($s * $s * $x2_sum) / ($n * $x2_sum - $x_sum * $x_sum));
  my $s_b1 = sqrt(($s * $s * $n) / ($n * $x2_sum - $x_sum * $x_sum));
  return $b0, $b1, $yfit, $s_b0, $s_b1;
}
