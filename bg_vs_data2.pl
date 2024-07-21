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
use PDL::Fit::Gaussian;
use PDL::Fit::Polynomial;
use Lab;

use Getopt::Long;
my %default_opts = (
		    bgfile => $Lab::MERGEDBG,
		    dev => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'dev=s', 'bgfile=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my @meds = (60, 70, 80, 90, 100, 110);
my $cmp_anode = 'C-Ka';

my @anodes = @ARGV ? @ARGV : @{+(Lab::lines())[0]};

dev $opts{dev}, 3,2;

print STDERR "Reading background data\n";
my $bg_hists = Lab::bg_hists($opts{bgfile});
my $bg_med = Lab::hists_median($bg_hists);
my $bg_n = $bg_hists->sumover;

my (%src_hists, %src_med, %src_n);

for my $anode ($cmp_anode, @anodes) {

  next if exists $src_hists{$anode};

  print STDERR "Reading $anode data\n";

  $src_hists{$anode} = Lab::src_hists($anode);
  $src_n{$anode} = $src_hists{$anode}->sumover;
  $src_med{$anode} = Lab::hists_median($src_hists{$anode});
#  $src_med{$anode} = +(Lab::hists_stats($src_hists{$anode}))[2];

#  bin hist $src_med{$anode}->where($src_med{$anode}!=0), -.5, 255.5, 1;

}

# now we look at C-Ka, for a set of medians, find all subtaps with
# those medians, then plot and source median in those subtaps with
# each anode versus the bg median in those subtaps

for my $med (@meds) {

  for my $anode ($cmp_anode, @anodes) {
    next if $anode eq $cmp_anode;

    my $i = whichND(
		    rint($src_med{$cmp_anode}) == rint($med)
		    );

    #
    # histogram of the background medians in these subtaps
    #
    my ($bg_med_histx, $bg_med_histy)= hist $bg_med->indexND($i)->where($bg_med->indexND($i) != 0), -0.5, 255.5, 1;

    #
    # fit a gaussian to the bg median distribution, use result to determine
    # which range of bg medians to compare with illuminated medians
    #
    my ($bg_med_cen, undef, $bg_med_fwhm) = fitgauss1d($bg_med_histx, $bg_med_histy);
    my ($bg_med_lo, $bg_med_hi) = map { ceil($bg_med_cen-$_), floor($bg_med_cen+$_) }
      $bg_med_fwhm/2.354*sqrt(2)*erfi(0.75);

    my $bg_meds_to_use = $bg_med_histx->where(
					      ($bg_med_histx>=$bg_med_lo) &
					      ($bg_med_histx<=$bg_med_hi)
					     );

    my (@x, @y, @n);
    for my $this_bg_med ($bg_meds_to_use->list) {
      my $ii = which(
		     ($bg_med->indexND($i) == $this_bg_med) &
		     ($src_med{$anode}->indexND($i) != 0)
		     );

      # FIXME: we should use the actual histogram data to arrive at a median
      my $src_med_mean = +($src_med{$anode}->indexND($i)->index($ii)->stats($src_n{$anode}->indexND($i)->index($ii)))[0];

      push @x, $this_bg_med;
      push @y, $src_med_mean;
      push @n, $src_n{$anode}->indexND($i)->index($ii)->sum;
    }

    my ($x, $y, $n) = (pdl(\@x), pdl(\@y), pdl(\@n));

  points($x, $y,
	 {
	  xtitle => 'background median',
	  ytitle => 'source data median',
	  title => 'C-Ka median = '.$med,
	 }
	);
  }
}

=begin comment

#bin $bg_med_histx, $bg_med_histy;

# hold fitted coefficients vs energy for each anode
my (@m, @b);

my @e = Lab::energies(@anodes);

for (0..$#anodes) {

  my $energy = $e[$_];
  my $anode = $anodes[$_];

  my (@x, @y, @n);
  # read in source data
  $src_hists .= 0;
  my ($src_med_histx, $src_med_histy)= hist $src_med->where($src_med != 0), -0.5, 255.5, 1;

    # FIXME: we should use the actual histogram data to arrive at a median
    my $src_med_mean = +($src_med->indexND($i)->stats($src_n->indexND($i)))[0];
    push @x, $this_bg_med;
    push @y, $src_med_mean;
    push @n, $src_n->indexND($i)->sum;
  }

  my ($x, $y, $n) = (pdl(\@x), pdl(\@y), pdl(\@n));

  points($x, $y,
	 {
	  xtitle => 'background median',
	  ytitle => 'source data median',
	  title => $anode,
	 }
	);
  hold;

  my ($yfit, $coeffs) = fitpoly1d($x, $y, 2, { Weights => $n });
  line $x, $yfit;
  hold;

  release;

  push @b, $coeffs->at(0);
  push @m, $coeffs->at(1);

}

$_/=1e3 for @e;

points pdl(\@e)->log10, pdl(\@m), { xtitle => 'energy', ytitle => 'slope', axis=>'logx' };
points pdl(\@e)->log10, pdl(\@b), { xtitle => 'energy', ytitle => 'intercept', axis=>'logx' };

=cut


exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
