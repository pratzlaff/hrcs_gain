#! /usr/bin/perl -w

use PDL;

use Inline Pdlpp;

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
use Chandra::Tools::Common qw( read_bintbl_cols );
use IO::Handle;
use File::Path qw( mkpath );
use Data::Dumper;
use Math::Trig qw( pi );
use PDL::Graphics::PGPLOT;
use PDL::Fit::LM qw( lmfit );
use Lab;

use Getopt::Long;

my $rdbext = '.rdb';
my $binext = '.bin';

my %default_opts = (
		    filtdir => $Lab::ANALDIR,
		    subtaps => 3,
		    outdir => '.',
		    rdb => 1,
		    bin => 1,
		    fit => 0,
		    fitcnts => 50, # subtap count threshold for fitting
		    bincnts => 1,  #   "      "       "      "  binfile output
		    trim => 5, # trimmed mean level
		    config => 'hrcs_lab.rdb',
		    twogauss => 0, # fit two normal distributions
		    subext => 1,
		    bgsubtract => 0,
		    bgfulltap => 0,
		    bgfile => $Lab::ANALDIR.'/merged_bg_evt1_filt_spi.fits',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'htmldir=s', 'filtdir=s', 'subtaps=i', 'subext!',
	   'binext=s', 'rdbext=s', 'rdb!', 'bin!', 'config=s',
	   'outdir=s', 'twogauss!', 'fit!',
	   'rdbfile=s', 'phabinfile=s', 'sampbinfile=s',
	   'spimeanbinfile=s', 'spimedbinfile=s',
	   'bgsubtract!',
	   'bgfile=s', 'bgfulltap!', '3x3!', 'bincnts=i',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

$opts{trim} /= 100;

if ($opts{'3x3'}) {
  $rdbext = '_3x3' . $rdbext;
  $binext  = '_3x3' . $binext;
}
$rdbext = $opts{rdbext} if $opts{rdbext};
$binext = $opts{binext} if $opts{binext};

#print join(', ', map(subtap_offsets($_), 0..2)),"\n";
#print join(', ', map(subtap_offsets2($_), 0..2)),"\n";
#exit;

$Lab::TESTFILE = $opts{config};
#my ($lines, undef, undef, undef, $hrc_files) = Lab::test_data() or die;
my ($lines, $hrc_files) = (Lab::test_data())[0,4] or die;

# if given basename arguments, use them, otherwise process all files
my @base = @ARGV;
if (!@base) {
  @base = @$hrc_files;
}

my @tmp;
my %lines;
@lines{@$lines} = ();
for my $b (@base) {
  if (exists $lines{$b}) {
    push @tmp, @{$hrc_files}[grep $lines->[$_] eq $b, 0..$#{$lines}];
  }
  else {
    push @tmp, $b;
  }
}
@base = @tmp;
for my $b (@base) {
  my $evt1 = "$opts{filtdir}/${b}_evt1_filt_spi.fits";
  -f $evt1 or die "could not find file '$evt1'";
}

my $exptime;
$exptime = (Lab::test_data_hrcfile(@base))[3] if $opts{bgsubtract};

my %bgtimes = Lab::merged_bg_exptimes();

#
# read merged background file
#
my @bg_cols = qw( rawx rawy pha samp spimean spimed chipid );
my %bg;
$bg{$_} = float([]) for @bg_cols;
if ($opts{bgsubtract}) {
  -f $opts{bgfile} or die "could not find file '$opts{bgfile}'";
  print STDERR "reading $opts{bgfile}...";
  STDERR->flush;
  (@bg{qw( rawx rawy pha samp spimean spimed chipid )}) = read_bintbl_cols($opts{bgfile}, qw( rawx rawy pha samp spimean spimed chip_id ), { extname => 'events', status => 1 } ) or die;
  print STDERR " done\n";
}
my %bg_wide = %bg;

my $subtaps = $opts{subtaps};

=begin comment

my ($xtap, $xsubtap, $ytap, $ysubtap, $x1, $x2, $y1, $y2) =
  mkcoords(Lab::RAWX_MIN, Lab::RAWX_MAX, Lab::RAWY_MIN, Lab::RAWY_MAX, Lab::TAPSIZE, $opts{subtaps});

=cut

#print $_->nelem,"\n" for $ytap, $ysubtap, $xtap, $xsubtap;
#wcols $ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2;
#exit;

my $xtap_reduced = sequence(long, long(Lab::RAWX_MAX)/long(Lab::TAPSIZE));
my $xsubtap_reduced = sequence(long, $subtaps);

my $ytap_reduced = sequence(long, long(Lab::RAWY_MAX)/long(Lab::TAPSIZE));
my $ysubtap_reduced = sequence(long, $subtaps);

# my $ytap_reduced = long(112);
# my $ysubtap_reduced = long(2);
# my $xtap_reduced = long(8,9);
# my $xsubtap_reduced = long(0,1,2);

for my $i (0..$#base) {
  my $base = $base[$i];
  my $exptime = $exptime->[$i];

  my $rdbfile = $opts{rdbfile} || "$opts{outdir}/$base$rdbext";
  my $pha_binfile = $opts{phabinfile} || "$opts{outdir}/${base}_pha$binext";
  my $samp_binfile = $opts{sampbinfile} || "$opts{outdir}/${base}_samp$binext";
  my $spimean_binfile = $opts{spimeanbinfile} || "$opts{outdir}/${base}_spimean$binext";
  my $spimed_binfile = $opts{spimedbinfile} || "$opts{outdir}/${base}_spimed$binext";

  my $rdb;
  if ($opts{rdb}) {
    $rdb = IO::File->new("> $rdbfile") or die "could not open $rdbfile: $!";
    my ($cols, $types) = rdb_colnames();
    $rdb->print( join("\t",@$cols), "\n", join("\t", @$types), "\n" );
  }

  my ($pha_bin, $samp_bin, $spimean_bin, $spimed_bin);
  if ($opts{bin}) {
    $pha_bin = Lab::OutBinFile->new($pha_binfile, 'pha') or die;
    $samp_bin = Lab::OutBinFile->new($samp_binfile, 'samp') or die;
    $spimean_bin = Lab::OutBinFile->new($spimean_binfile, 'spimean') or die;
    $spimed_bin = Lab::OutBinFile->new($spimed_binfile, 'spimed') or die;
  }

  my @src_cols = qw( rawx rawy pha samp spimean spimed );
  my %src;

  my $evt1 = "$opts{filtdir}/${base}_evt1_filt_spi.fits";
  -f $evt1 or die $base;

  print STDERR "reading $evt1...";
  STDERR->flush;
  (@src{qw( rawx rawy pha samp spimean spimed )}) = read_bintbl_cols($evt1, 'rawx', 'rawy', 'pha', 'samp', 'spimean', 'spimed', { extname => 'events', status => 1 } ) or die;
  print STDERR " done\n";

  my %src_wide = %src;

=begin comment

  for my $i (0..$xtap->nelem-1) {
    my ($xtap, $xsubtap, $ytap, $ysubtap, $x1, $x2, $y1, $y2) =
      (
       $xtap->at($i),
       $xsubtap->at($i),
       $ytap->at($i),
       $ysubtap->at($i),
       $x1->at($i),
       $x2->at($i),
       $y1->at($i),
       $y2->at($i),
       );

    my $index = which(
		      ($rawx >= $x1) & ($rawx <= $x2) &
		      ($rawy >= $y1) & ($rawy <= $y2)
		      );
    $index->nelem or next;

  }

=cut


  # tried filtering at various steps and changing the order of rawx/y
  # stepping, the following represents my best effort

  for my $ytap ($ytap_reduced->list) {
    my %bg = %bg;
    my %bg_wide = %bg_wide;
    my %src = %src;
    my %src_wide = %src_wide;

    my ($y1, $y2) = tap_raw_range($ytap);

    # save enough data for inclusion of subtaps on either side
    my $yw1 = (subtap_raw_range($ytap-1, $subtaps-1))[0];
    my $yw2 = (subtap_raw_range($ytap+1, 0))[1];

    my $i = which( ($src_wide{rawy} >= $yw1) & ($src_wide{rawy} <= $yw2) );
    $src_wide{$_} = $src_wide{$_}->index($i) for @src_cols;

    $i = which( ($src{rawy} >= $y1) & ($src{rawy} <= $y2) );
    $src{$_} = $src{$_}->index($i) for @src_cols;


    if ($opts{bgsubtract}) {
      if (!$opts{bgfulltap}) {
	$i = which( ($bg_wide{rawy} >= $yw1) & ($bg_wide{rawy} <= $yw2) );
	$bg_wide{$_} = $bg_wide{$_}->index($i) for @bg_cols;
      }

      $i = which( ($bg{rawy} >= $y1) & ($bg{rawy} <= $y2) );
      $bg{$_} = $bg{$_}->index($i) for @bg_cols;
    }

    for my $ysubtap ($ysubtap_reduced->list) {
      my %bg = %bg;
      my %bg_wide = %bg_wide;
      my %src = %src;
      my %src_wide = %src_wide;

      my ($y1, $y2) = subtap_raw_range($ytap, $ysubtap);

      my $i = which( ($src{rawy} >= $y1) & ($src{rawy} <= $y2) );
      $src{$_} = $src{$_}->index($i) for @src_cols;

      if ($opts{bgsubtract} and !$opts{bgfulltap}) {
	$i = which( ($bg{rawy} >= $y1) & ($bg{rawy} <= $y2) );
	$bg{$_} = $bg{$_}->index($i) for @bg_cols;
      }

      for my $xtap ($xtap_reduced->list) {
	my %bg = %bg;
	my %bg_wide = %bg_wide;
	my %src = %src;
	my %src_wide = %src_wide;

	my ($x1, $x2) = tap_raw_range($xtap);

	# save enough data for inclusion of subtaps on either side
	my $xw1 = (subtap_raw_range($xtap-1, $subtaps-1))[0];
	my $xw2 = (subtap_raw_range($xtap+1, 0))[1];

	my $i = which( ($src_wide{rawx} >= $xw1) & ($src_wide{rawx} <= $xw2) );
	$src_wide{$_} = $src_wide{$_}->index($i) for @src_cols;

	$i = which( ($src{rawx} >= $x1) & ($src{rawx} <= $x2) );
	$src{$_} = $src{$_}->index($i) for @src_cols;

	if ($opts{bgsubtract}) {

	  if (!$opts{bgfulltap}) {
	    $i = which( ($bg_wide{rawx} >= $xw1) & ($bg_wide{rawx} <= $xw2) );
	    $bg_wide{$_} = $bg_wide{$_}->index($i) for @bg_cols;
	  }

	  $i = which( ($bg{rawx} >= $x1) & ($bg{rawx} <= $x2) );
	  $bg{$_} = $bg{$_}->index($i) for @bg_cols;
	}

	for my $xsubtap ($xsubtap_reduced->list) {
	  my %bg = %bg;
	  my %bg_wide = %bg_wide;
	  my %src = %src;
	  my %src_wide = %src_wide;

	  my ($x1, $x2) = subtap_raw_range($xtap, $xsubtap);

	  # might have to modify these for low counts
	  my ($y1, $y2) = ($y1->copy, $y2->copy);

	  my $i = which( ($src{rawx} >= $x1) & ($src{rawx} <= $x2) );
	  my $norig = $i->nelem;
 	  my $spimean = $src{spimean}->index($i);
 	  my $spimed = $src{spimed}->index($i);
 	  my $samp = $src{samp}->index($i);
 	  my $pha = $src{pha}->index($i);

	  my ($spimean_bg, $spimed_bg, $samp_bg, $pha_bg, $chipid_bg) = @bg{qw( spimean spimed samp pha chipid )};
	  if ($opts{bgsubtract} and !$opts{bgfulltap}) {
	    my $i = which( ($bg{rawx} >= $x1) & ($bg{rawx} <= $x2) );
	    $spimean_bg = $bg{spimean}->index($i);
	    $spimed_bg = $bg{spimed}->index($i);
	    $samp_bg = $bg{samp}->index($i);
	    $pha_bg = $bg{pha}->index($i);
	    $chipid_bg = $bg{chipid}->index($i);
	  }

	  my $subs = '1x1';

	  if ($opts{subext} or $opts{'3x3'}) {

	    if ($norig < 100 or $opts{'3x3'}) {
	      $subs = '3x3';

	      $y1 = (subtap_raw_range($ytap, $ysubtap-1))[0];
	      $y2 = (subtap_raw_range($ytap, $ysubtap+1))[1];
	      $x1 = (subtap_raw_range($xtap, $xsubtap-1))[0];
	      $x2 = (subtap_raw_range($xtap, $xsubtap+1))[1];

# 	      $y1 -= long(Lab::TAPSIZE/$subtaps);
# 	      $y2 = $y1 + long(Lab::TAPSIZE * (3 /$subtaps)) - 1;
# 	      $x1 -= long(Lab::TAPSIZE/$subtaps);
# 	      $x2 = $x1 + long(Lab::TAPSIZE * (3 /$subtaps)) - 1;

 	      $x1 = Lab::RAWX_MIN if $x1 < Lab::RAWX_MIN;
 	      $x2 = Lab::RAWX_MAX if $x2 > Lab::RAWX_MAX;
 	      $y1 = Lab::RAWY_MIN if $y1 < Lab::RAWY_MIN;
 	      $y2 = Lab::RAWY_MAX if $y2 > Lab::RAWY_MAX;

	      my $i = which(
			    ($src_wide{rawx} >= $x1) &
			    ($src_wide{rawx} <= $x2) &
			    ($src_wide{rawy} >= $y1) &
			    ($src_wide{rawy} <= $y2)
			   );
	      $spimean = $src_wide{spimean}->index($i);
	      $spimed = $src_wide{spimed}->index($i);
	      $samp = $src_wide{samp}->index($i);
	      $pha = $src_wide{pha}->index($i);

	      if ($opts{bgsubtract} and !$opts{bgfulltap}) {
		$i = which(
			   ($bg_wide{rawx} >= $x1) &
			   ($bg_wide{rawx} <= $x2) &
			   ($bg_wide{rawy} >= $y1) &
			   ($bg_wide{rawy} <= $y2)
			  );
		$spimean_bg = $bg_wide{spimean}->index($i);
		$spimed_bg = $bg_wide{spimed}->index($i);
		$samp_bg = $bg_wide{samp}->index($i);
		$pha_bg = $bg_wide{pha}->index($i);
		$chipid_bg = $bg_wide{chipid}->index($i);
	      }

	    } elsif ($norig < 150 and !$opts{'3x3'}) {
	      $subs = '3x1';


	      $y1 = (subtap_raw_range($ytap, $ysubtap-1))[0];
	      $y2 = (subtap_raw_range($ytap, $ysubtap+1))[1];

#  	      $y1 -= long(Lab::TAPSIZE/$subtaps);
#  	      $y2 = $y1 + long(Lab::TAPSIZE * (3 /$subtaps)) - 1;

 	      $y1 = Lab::RAWY_MIN if $y1 < Lab::RAWY_MIN;
 	      $y2 = Lab::RAWY_MAX if $y2 > Lab::RAWY_MAX;

	      my $i = which(
			    ($src_wide{rawx} >= $x1) &
			    ($src_wide{rawx} <= $x2) &
			    ($src_wide{rawy} >= $y1) &
			    ($src_wide{rawy} <= $y2)
			   );
	      $spimean = $src_wide{spimean}->index($i);
	      $spimed = $src_wide{spimed}->index($i);
	      $samp = $src_wide{samp}->index($i);
	      $pha = $src_wide{pha}->index($i);

	      if ($opts{bgsubtract} and !$opts{bgfulltap}) {
		$i = which(
			   ($bg_wide{rawx} >= $x1) &
			   ($bg_wide{rawx} <= $x2) &
			   ($bg_wide{rawy} >= $y1) &
			   ($bg_wide{rawy} <= $y2)
			  );
		$spimean_bg = $bg_wide{spimean}->index($i);
		$spimed_bg = $bg_wide{spimed}->index($i);
		$samp_bg = $bg_wide{samp}->index($i);
		$pha_bg = $bg_wide{pha}->index($i);
		$chipid_bg = $bg_wide{chipid}->index($i);
	      }

	    }
	  }

	  my $n = $pha->nelem;
	  my $nnet = $n;

	  die $samp->nelem . ' ' . $pha->nelem
	    unless $samp->nelem == $pha->nelem;

	  die $samp->nelem . ' ' . $spimean->nelem
	    unless $samp->nelem == $spimean->nelem;

	  die $samp->nelem . ' ' . $spimed->nelem
	    unless $samp->nelem == $spimed->nelem;

	  my ($pha_mean, $pha_tmean, $pha_rms, $pha_median, $pha_iqr,
	      $samp_mean, $samp_tmean, $samp_rms, $samp_median, $samp_iqr,
	      $spimean_mean, $spimean_tmean, $spimean_rms, $spimean_median, $spimean_iqr,
	      $spimed_mean, $spimed_tmean, $spimed_rms, $spimed_median, $spimed_iqr,
	     ) = (0)x20;

	  goto PRINTRDB unless $n;

	  if ($opts{bgsubtract} && $samp_bg->nelem) {

	    my ($pha_x, $pha_spec) = hist($pha->double, -0.5, 255.5, 1);
	    my ($samp_x, $samp_spec) = hist($samp, -0.5, 511.5, 1);
	    my ($spimean_x, $spimean_spec) = hist($spimean, -0.5, 511.5, 1);
	    my ($spimed_x, $spimed_spec) = hist($spimed, -0.5, 511.5, 1);

	    my $pha_bg_spec = zeroes($pha_spec->nelem);
	    my $samp_bg_spec = zeroes($samp_spec->nelem);
	    my $spimean_bg_spec = zeroes($spimean_spec->nelem);
	    my $spimed_bg_spec = zeroes($spimed_spec->nelem);

	    if ($opts{bgfulltap}) {
	      # FIXME: get the real chipid
	      my $w = $exptime / $bgtimes{1} * ($x2-$x1+1) * ($y2-$y1+1) / Lab::TAPSIZE / Lab::TAPSIZE;
	      $pha_bg_spec += $w * hist($pha_bg, -0.5, 255.5, 1);
	      $samp_bg_spec += $w * hist($samp_bg, -0.5, 511.5, 1);
	      $spimean_bg_spec += $w * hist($spimean_bg, -0.5, 511.5, 1);
	      $spimed_bg_spec += $w * hist($spimed_bg, -0.5, 511.5, 1);
	      $nnet -= $w * $pha_bg->nelem;
	    }
	    else {
	      for my $chipid ($chipid_bg->uniq->list) {
		my $i = which($chipid_bg == $chipid);
		my $w = $exptime / $bgtimes{$chipid};
		$pha_bg_spec += $w * hist($pha_bg->index($i), -0.5, 255.5, 1);
		$samp_bg_spec += $w * hist($samp_bg->index($i), -0.5, 511.5, 1);
		$spimean_bg_spec += $w * hist($spimean_bg->index($i), -0.5, 511.5, 1);
		$spimed_bg_spec += $w * hist($spimed_bg->index($i), -0.5, 511.5, 1);
		$nnet -= $w * $i->nelem;
	      }
	    }

	    $pha_spec -= $pha_bg_spec;
	    $samp_spec -= $samp_bg_spec;
	    $spimean_spec -= $spimean_bg_spec;
	    $spimed_spec -= $spimed_bg_spec;

	    ($pha_mean, $pha_tmean, $pha_rms, $pha_iqr, $pha_median) = hist_stats($pha_x, $pha_spec);
	    ($samp_mean, $samp_tmean, $samp_rms, $samp_iqr, $samp_median) = hist_stats($samp_x, $samp_spec);
	    ($spimean_mean, $spimean_tmean, $spimean_rms, $spimean_iqr, $spimean_median) = hist_stats($spimean_x, $spimean_spec);
	    ($spimed_mean, $spimed_tmean, $spimed_rms, $spimed_iqr, $spimed_median) = hist_stats($spimed_x, $spimed_spec);
	  }

	  # simple statistics if no background subtraction
	  else {
	    my $pha = $pha->qsort->double;
	    my $samp = $samp->qsort;
	    my $spimean = $spimean->qsort;
	    my $spimed = $spimed->qsort;
	    ($pha_mean, $pha_rms, $pha_median) = $pha->stats;
	    ($samp_mean, $samp_rms, $samp_median) = $samp->stats;
	    ($spimean_mean, $spimean_rms, $spimean_median) = $spimean->stats;
	    ($spimed_mean, $spimed_rms, $spimed_median) = $spimed->stats;
	    $spimean_tmean = Lab::trim($spimean, $opts{trim}, 1)->avg;
	    $spimed_tmean = Lab::trim($spimed, $opts{trim}, 1)->avg;
	    $samp_tmean = Lab::trim($samp, $opts{trim}, 1)->avg;
	    $pha_tmean = Lab::trim($pha, $opts{trim}, 1)->avg;
	    if ($pha->nelem > 3) {
	      $spimean_iqr = Lab::iqr($spimean, 1);
	      $spimed_iqr = Lab::iqr($spimed, 1);
	      $samp_iqr = Lab::iqr($samp, 1);
	      $pha_iqr = Lab::iqr($pha, 1);
	    }
	  }

	  my $pha_params = zeroes($opts{twogauss} ? 6 : 3);
	  my $samp_params = $pha_params->copy;
	  my $spimean_params = $pha_params->copy;
	  my $spimed_params = $pha_params->copy;

	  if ($pha->nelem >= $opts{bincnts} and $opts{bin}) {
	    $pha_bin->add_subtap($ytap, $ysubtap, $xtap, $xsubtap,
				 $y1, $y2, $x1, $x2, $pha);
	    $samp_bin->add_subtap($ytap, $ysubtap, $xtap, $xsubtap,
				  $y1, $y2, $x1, $x2, $samp);
	    $spimean_bin->add_subtap($ytap, $ysubtap, $xtap, $xsubtap,
				  $y1, $y2, $x1, $x2, $spimean);
	    $spimed_bin->add_subtap($ytap, $ysubtap, $xtap, $xsubtap,
				  $y1, $y2, $x1, $x2, $spimed);
	  }

	  if ($pha->nelem >= $opts{fitcnts} and $opts{fit}) {

	    my ($pha_x, $pha_y) = hist($pha->double, -0.5, 255.5, 1);
	    my ($samp_x, $samp_y) = hist($samp, -0.5, 511.5, 1);
	    my ($spimean_x, $spimean_y) = hist($spimean, -0.5, 511.5, 1);
	    my ($spimed_x, $spimed_y) = hist($spimed, -0.5, 511.5, 1);

	    my ($pha_fit, $samp_fit, $spimean_fit, $spimed_fit);

	    # params are norm, sigma, mean (x2)
	    my $pha_params_init = $opts{twogauss} ?
	      float($pha_y->sum/5, 50, $pha_median/2, $pha_y->sum, 15, $pha_median) :
		float($pha_y->sum, 15, $pha_median);
	    my $samp_params_init = $pha_params_init->copy;
	    my $spimean_params_init = $pha_params_init->copy;
	    my $spimed_params_init = $pha_params_init->copy;

	    ($pha_fit, $pha_params) = lmfit($pha_x, $pha_y, 1, $opts{twogauss} ? \&two_normals : \&one_normal, $pha_params_init);

	    ($samp_fit, $samp_params) = lmfit($samp_x, $samp_y, 1, $opts{twogauss} ? \&two_normals : \&one_normal, $samp_params_init);

	    ($spimean_fit, $spimean_params) = lmfit($spimean_x, $spimean_y, 1, $opts{twogauss} ? \&two_normals : \&one_normal, $spimean_params_init);

	    ($spimed_fit, $spimed_params) = lmfit($spimed_x, $spimed_y, 1, $opts{twogauss} ? \&two_normals : \&one_normal, $spimed_params_init);

	  }

	PRINTRDB:
	  if ($opts{rdb}) {
	    $rdb->print(join("\t",
		     $ytap, $ysubtap,
		     $xtap, $xsubtap, "$y1:$y2", "$x1:$x2", $subs,
		     ($opts{bgsubtract} ?
		      (
		       sprintf("%.1f", $nnet),
		       $norig, $n,
		       sprintf("%.1f", $n-$nnet),
		      ) :
		      (
		       $n, $norig,
			which($pha<3)->nelem,
		       which($pha==255)->nelem
		      )
		     ),

#		     sprintf("%.1f", $pha_mean),
		     sprintf("%.1f", $pha_tmean),
		     sprintf("%.1f", $pha_rms),
		     ($opts{bgsubtract} ? sprintf("%.1f", $pha_median) : $pha_median),

#		     sprintf("%.1f", $samp_mean),
		     sprintf("%.1f", $samp_tmean),
		     sprintf("%.1f", $samp_rms),
		     sprintf("%.1f", $samp_median),

#		     sprintf("%.1f", $spimean_mean),
		     sprintf("%.1f", $spimean_tmean),
		     sprintf("%.1f", $spimean_rms),
		     sprintf("%.1f", $spimean_median),

#		     sprintf("%.1f", $spimed_mean),
		     sprintf("%.1f", $spimed_tmean),
		     sprintf("%.1f", $spimed_rms),
		     sprintf("%.1f", $spimed_median),

		     ( $opts{fit} ? (
		       map { sprintf "%.2f", $_ } $pha_params->list,
		       map { sprintf "%.2f", $_ } $samp_params->list,
		       map { sprintf "%.2f", $_ } $spimean_params->list,
		       map { sprintf "%.2f", $_ } $spimed_params->list,
				    ) : ()
		     ),

		     sprintf("%.1f", $pha_iqr),
		     sprintf("%.1f", $samp_iqr),
		     sprintf("%.1f", $spimean_iqr),
		     sprintf("%.1f", $spimed_iqr),
		    ), "\n");
	  }

	} # xsubtap
      }   # xtap
    }     # ysubtap
  }       # ytap



}

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

sub mkcoords {
  my ($xmin, $xmax, $ymin, $ymax, $tapsize, $subtaps) = @_;

  my $xdim = $xmax - $xmin + 1;
  my $ydim = $ymax - $ymin + 1;

  my $xtap_reduced = sequence(long, long($xdim)/long($tapsize));
  my $xsubtap_reduced = sequence(long, $subtaps);

  my $ytap_reduced = sequence(long, long($ydim)/long($tapsize));
  my $ysubtap_reduced = sequence(long, $subtaps);

  # expand rawX tap/subtap lists to span a single rawY subtap
  my $xtap = $xtap_reduced->dummy(0, $xsubtap_reduced->nelem)->flat;
  my $xsubtap = $xsubtap_reduced->dummy(-1, $xtap_reduced->nelem)->flat;

  # expand rawY tap/subtap lists fully
  my $ytap = $ytap_reduced->dummy(0, $xtap->nelem * $ysubtap_reduced->nelem)->flat;
  my $ysubtap = $ysubtap_reduced->dummy(0, $xtap->nelem)->flat->dummy(-1, $ytap_reduced->nelem)->flat;

  # finish expansion of rawX tap/subtap lists
  $_ = $_->dummy(-1, $ytap_reduced->nelem * $ysubtap_reduced->nelem)->flat
    for $xtap, $xsubtap;

  my $x1 = $xtap*$tapsize+Lab::RAWX_MIN + long($tapsize*$xsubtap/$subtaps);
  my $x2 = $xtap*$tapsize+Lab::RAWX_MIN + long($tapsize*($xsubtap+1)/$subtaps) - 1;

  my $y1 = $ytap*$tapsize+Lab::RAWY_MIN + long($tapsize*$ysubtap/$subtaps);
  my $y2 = $ytap*$tapsize+Lab::RAWY_MIN + long($tapsize*($ysubtap+1)/$subtaps) - 1;

  return $xtap, $xsubtap, $ytap, $ysubtap, $x1, $x2, $y1, $y2;
}

sub hist_stats {
  my ($x, $y) = @_;

  my $y_sum = $y->sum;

  return (0)x5 unless $y->sum > 0;

  my $mean = ($x*$y)->sum / $y_sum;

  my $diff = $x-$mean;

  my $y_copy = $y->copy;
  $y_copy->shift_negs;
  my $i = which ($y_copy >= 0);
  my $rms = $y_copy->index($i)->sum <=1 ?
    0 : sqrt(($y_copy->index($i) *  $diff->index($i) * $diff->index($i))->sum / ($y_copy->index($i)->sum-1));

  my ($median, $lq, $uq) = (0)x3;
  eval {
    my $cumu_bounds = zeroes($x->nelem);
     (my $tmp = $cumu_bounds->slice('0:-2')) .= ($x->slice('0:-2') + $x->slice('1:-1'))/2;
    $cumu_bounds->set(-1, $x->at(-1)+($x->at(-1)-$x->at(-2))/2);
    my $cumu = $y->cumusumover / $y_sum;
    $uq = interpol(0.75, $cumu, $cumu_bounds);
    $median = interpol(0.5, $cumu, $cumu_bounds);
    $lq = interpol(0.25, $cumu, $cumu_bounds);
  };

  my $trimmed_y = Lab::hists_trim($y_copy, $opts{trim});
  my $trimmed_y_sum = $trimmed_y->sum;
  my $tmean = $trimmed_y_sum > 0 ? ($x*$trimmed_y)->sum / $trimmed_y->sum : 0;

  return $mean, $tmean, $rms, $uq-$lq, $median;
}

sub subtap_offsets_old {
  my $subtap = shift;
  my $o1 = long(Lab::TAPSIZE*$subtap/$opts{subtaps});
  my $o2 = long(Lab::TAPSIZE*($subtap+1)/$opts{subtaps}) - 1 ;
  return $o1, $o2;
}

sub subtap_offsets {
  my $subtap_orig = shift;

  my $subtap = $subtap_orig;

  $subtap = $subtap + $opts{subtaps} if $subtap_orig < 0;

  my $o1 = rint(Lab::TAPSIZE*$subtap/$opts{subtaps});
  my $o2 = rint(Lab::TAPSIZE*($subtap+1)/$opts{subtaps}) - 1 ;

  $o1-=Lab::TAPSIZE, $o2-=Lab::TAPSIZE if $subtap_orig < 0;

  return $o1, $o2;
}

sub tap_raw_range {
  my $tap = shift;
  my $c1 = $tap * Lab::TAPSIZE + Lab::RAWX_MIN; # RAWX_MIN == RAWY_MIN
  my $c2 = $c1 + Lab::TAPSIZE - 1;
  return $c1, $c2;
}

sub subtap_raw_range {
  my ($tap, $subtap) = @_;
  my $tc1 = (tap_raw_range($tap))[0]; # begin coord of tap
  my ($stc1, $stc2) = subtap_offsets($subtap);
  return $tc1+$stc1, $tc1+$stc2;
}

sub rdb_colnames {
  my @cols = qw( crsv vsub crsu usub rawy_range rawx_range subs);
  push @cols, ( $opts{bgsubtract} ?
		qw( nnet norig n BG ) :
		qw (  n norig pha_lt_3 pha255 )
	      );
  my $trim = $opts{trim} * 100;
  push @cols,
#    'pmean',
    'ptmean'.$trim, qw( prms pmed ),
#    'smean',
    'stmean'.$trim, qw( srms smed ),
#    'spimeanmean',
    'spimeantmean'.$trim, qw( spimeanrms spimeanmed ),
#    'spimedmean',
    'spimedtmean'.$trim, qw( spimedrms spimedmed ),
     qw( piqr siqr spimeaniqr spimediqr );

  my @fit_cols = qw( _gnorm _gsigma _gmean );

  if ($opts{fit}) {
    push @cols, map('pha'.$_, @fit_cols);
    push @cols, map('pha'.$_.'2', @fit_cols) if $opts{twogauss};

    push @cols, map('samp'.$_, @fit_cols);
    push @cols, map('samp'.$_.'2', @fit_cols) if $opts{twogauss};

    push @cols, map('spimean'.$_, @fit_cols);
    push @cols, map('spimean'.$_.'2', @fit_cols) if $opts{twogauss};

    push @cols, map('spimed'.$_, @fit_cols);
    push @cols, map('spimed'.$_.'2', @fit_cols) if $opts{twogauss};
  }

  my @types = ( qw( N N N N S S S ), ('N')x(@cols-7) );

  return \(@cols, @types);
}

__DATA__

__Pdlpp__

pp_def('shift_negs',
	 Pars => 'i(n);',
       Code => 'int i, n_size; n_size = $SIZE(n); for (i=0; i<n_size-1; ++i) { if ($i(n=>i) < 0) { $i(n=>i+1) += $i(n=>i); $i(n=>i) = 0; } }',
	);
