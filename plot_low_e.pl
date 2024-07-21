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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> April 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use PGPLOT;
use MyRDB;
use Data::Dumper;
use Chandra::Tools::Common qw();

use Getopt::Long;
my %default_opts = (
		    rdb => 'hrcs_lab.rdb',
		    datadir => '/data/legs/rpete/data/hrcs_lab/analysis',
		    mincnts => 50,
		    dev => '/xs',
		    charsize => 2,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'rdb=s', 'datadir=s', 'mincnts=i', 'dev=s',
	   'revcolors!', 'colors!', 'charsize=f',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my ($line, $energy, $hrc_file) = MyRDB::rdb_cols($opts{rdb}, qw( line energy HRC_file )) or die;

my @anodes = @ARGV ? @ARGV : @$line;

# remove duplicates
my %anodes;
@anodes{@anodes} = ();
@anodes = keys %anodes;

# ensure each exists in our rdb file
my %lines;
@lines{@$line} = @$energy;

for (@anodes) {
  exists $lines{$_} or die "no entries in $opts{rdb} for anode '$_'";
}

# order by energy
@anodes = sort { $lines{$a} <=> $lines{$b} } @anodes;

dev $opts{dev}, 3, 3;

# y limits on end of each plate, 64 taps per plate, 256 pixels per tap
my @ylim_hi = map { 64 * 256 * $_ } 1..3;
my @ylim_lo = map { $_ - 64 * 256 + 1 } @ylim_hi;
my @colors = (1,2,3);

if ($opts{revcolors}) {
  @colors = reverse @colors;
  @ylim_hi = reverse @ylim_hi;
  @ylim_lo = reverse @ylim_lo;
}

# median/std dev, in each energy, for various quantities of interest
my %stats = ();
# example entry
#    'mupeak/PHAmedian' => { energy => [], mean => [], prms => [] },


for my $anode (@anodes) {

  my @i = grep { $line->[$_] eq $anode } 0..$#{$line};

  # energy (eV) for this anode
  my $energy = $energy->[$i[0]];

  my @hrc_file = @{$hrc_file}[@i];

  my $n = long [];
  my $y1 = long [];
  my $y2 = $y1->copy;
  my $mean = float [];
  my $median = $mean->copy;
  my $gnorm = $mean->copy;
  my $gsigma = $mean->copy;
  my $gmean = $mean->copy;
  my $resmean = $mean->copy;
  my $resrms = $mean->copy;

  for my $hrc_file (@hrc_file) {
    my $rdb = $opts{datadir} . '/' . $hrc_file . '.rdb';
    my $rdb_lowe = $opts{datadir} . '/' . $hrc_file . '_lowe.rdb';
    for ($rdb, $rdb_lowe) {
      -f $_ or die "'$_' does not exist";
    }

    my ($tn, $tmean, $tmedian)
      = MyRDB::rdb_cols( $rdb, qw( n mean median ) )
	or die;

    my ($ty1, $ty2, $tsum, $tgnorm, $tgsigma, $tgmean, $tresmean, $tresrms)
      = MyRDB::rdb_cols( $rdb_lowe, qw( y1 y2 sum gnorm gsigma gmean resmean resrms ) )
	or die;

    $_ = long $_ for $tn, $ty1, $ty2, $tsum;
    $_ = float $_ for $tmean, $tmedian, $tgnorm, $tgsigma, $tgmean, $tresmean, $tresrms;

    my $i1 = which($tn >= $opts{mincnts});
    my $i2 = which($tsum >= $opts{mincnts});

    next unless $i1->nelem;

    $i1->nelem == $i2->nelem or die $hrc_file;

    # ensure we're looking at the same subtaps in both files
    my $sum_diff = $tn->index($i1) - $tsum->index($i2);
    $sum_diff->abs->max == 0 or die $hrc_file;

    $n = $n->append($tn->index($i1));
    $mean = $mean->append($tmean->index($i1));
    $median = $median->append($tmedian->index($i1));

    $y1 = $y1->append($ty1->index($i2));
    $y2 = $y2->append($ty2->index($i2));
    $gnorm = $gnorm->append($tgnorm->index($i2));
    $gsigma = $gsigma->append($tgsigma->index($i2));
    $gmean = $gmean->append($tgmean->index($i2));
    $resmean = $resmean->append($tresmean->index($i2));
    $resrms = $resrms->append($tresrms->index($i2));
  }

  my %items = (
	       'mupeak/PHAmedian' => $gmean / $median,
	       'PHAmedian/PHAmean' => $median / $mean,
	       'mupeak/PHAmean' => $gmean / $mean,

	       'mures/PHAmean' => $resmean / $mean,
	       'mures/PHAmedian' => $resmean / $median,
	       'mures/mupeak' => $resmean / $gmean,

	       'sigmapeak/PHAmean' => $gsigma / $mean,
	       'sigmapeak/PHAmedian' => $gsigma / $median,
	       'sigmapeak/mupeak' => $gsigma / $gmean,

	       'sigmares/PHAmean' => $resrms / $mean,
	       'sigmares/PHAmedian' => $resrms / $median,
	       'sigmares/mupeak' => $resrms / $gmean,

	       'Nres/Ntotal' => +($n-$gnorm)/$n,
	      );

  for my $key (keys %items) {
    my ($mean, $prms) = $items{$key}->stats;
    push @{$stats{$key}{energy}}, $energy;
    push @{$stats{$key}{mean}}, $mean;
    push @{$stats{$key}{prms}}, $prms;
  }

  points $mean, $gmean/$median, { xtitle => 'PHA mean', ytitle => '\\gm\\dpeak\\u / PHA median', title => $anode, xrange => [50, 200], yrange => [0.8, 1.25], symbol => -1, charsize => $opts{charsize} };

  points $mean, $median/$mean, { xtitle => 'PHA mean', ytitle => 'PHA median / PHA mean', title => $anode, xrange => [50, 200], yrange => [0.80, 1.25], symbol => -1, charsize => $opts{charsize} };
  if ($opts{colors}) {
    hold;
    for (0..$#ylim_hi) {
      my $i = which( ( $y1 >= $ylim_lo[$_] ) &
		     ( $y2 <= $ylim_hi[$_] )
		   );
      $i->nelem or next;
      points $mean->index($i), $median->index($i)/$mean->index($i), { symbol => -1, color => $colors[$_] };
    }
    release;
  }

  points $mean, $gmean/$mean, { xtitle => 'PHA mean', ytitle => '\\gm\\dpeak\\u / PHA mean', title => $anode, xrange => [50, 200], yrange => [0.80, 1.25], symbol => -1, charsize => $opts{charsize} };

  points $mean, $resmean/$mean, { xtitle => 'PHA mean', ytitle => '\\gm\\dres\\u / PHA mean', title => $anode, xrange => [50, 200], yrange => [0.30, 0.80], symbol => -1, charsize => $opts{charsize} };

  points $n, $resmean/$gmean, { ytitle => '\\gm\\dres\\u / \\gm\\dpeak\\u', xtitle => 'Total counts', title => $anode, xrange => [0, 1500], yrange => [0, 0.8], symbol => -1, charsize => $opts{charsize} };

  points $mean, $gsigma/$gmean, { xtitle => 'PHA mean', ytitle => '\\gs\\dpeak\\u / \\gm\\dpeak\\u', title => $anode, xrange => [50, 200], yrange => [0.05, 0.35], symbol => -1, charsize => $opts{charsize} };

  points $mean, $resrms/$resmean, { xtitle => 'PHA mean', ytitle => '\\gs\\dres\\u / \\gm\\dres\\u', title => $anode, xrange => [50, 200], yrange => [0.1, 0.6], symbol => -1, charsize => $opts{charsize} };

  points $n, +($n-$gnorm)/$n, { ytitle => 'Residual fractional counts', xtitle => 'Total counts', title => $anode, xrange => [0, 1500], yrange => [0, 0.4], symbol => -1, charsize => $opts{charsize} };

  points $mean, +($n-$gnorm)/$n, { ytitle => 'Residual fractional counts', xtitle => 'PHA mean', title => $anode, xrange => [50, 200], yrange => [0, 0.4], symbol => -1, charsize => $opts{charsize} };


=begin comment

  points $gsigma, $resrms, { xtitle => '\\gs\\dpeak', ytitle => '\\gs\\dres', symbol => -1 };
  bin hist $resrms; #, { xtitle => '\\gs\\dres', ytitle => 'n' };
  pglabel('\\gs\\dres', 'N', '');
#  points $g1pos, $g2pos, { xtitle => '\\gm\\d1\\u', ytitle => '\\gm\\d2\\u', title => $anode, xrange => [0,255], yrange => [0,255], symbol => -1 };

=cut

}

# convert all statistics to piddles
for my $key (keys %stats) {
  $stats{$key}{$_} = float $stats{$key}{$_} for qw( energy mean prms );
  $stats{$key}{energy} = $stats{$key}{energy}->log10 - 3;
}

my @keys = (
	    [ 'mupeak/PHAmedian', 'PHAmedian/PHAmean', 'mupeak/PHAmean' ],
	    [ 'mures/PHAmean', 'mures/PHAmedian', 'mures/mupeak' ],
	    [ 'sigmapeak/PHAmean', 'sigmapeak/PHAmedian', 'sigmapeak/mupeak' ],
	    [ 'sigmares/PHAmean', 'sigmares/PHAmedian', 'sigmares/mupeak' ],
	    [ 'Nres/Ntotal' ],
	    );

my @legends = (
	       ['\\gm\\dpeak\\u / PHA median', 'PHA median / PHA mean', '\\gm\\dpeak\\u / PHA mean'],
	       ['\\gm\\dres\\u / PHA mean', '\\gm\\dres\\u / PHA median', '\\gm\\dres\\u / \\gm\\dpeak\\u'],
	       ['\\gs\\dpeak\\u / PHA mean', '\\gs\\dpeak\\u / PHA median', '\\gs\\dpeak\\u / \\gm\\dpeak\\u'],
	       ['\\gs\\dres\\u / PHA mean', '\\gs\\dres\\u / PHA median', '\\gs\\dres\\u / \\gm\\dpeak\\u'],
	       ['N\\dres\\u / N\\dtotal\\u' ],
	       );

# plot of statistic vs energy for the various ratios
pgpage();
pgsubp(3,2);
for my $i (0..$#keys) {
  my @keys = @{$keys[$i]};
  my @legends = @{$legends[$i]};
  my @symbols = (0, 2, 3);
  my @colors = 1..3;

  # chop off unnecessary elements of arrays
  @$_ = @{$_}[0..$#keys] for \@symbols, \@colors;

  # get plot limits
  my ($elo, $ehi, $ylo, $yhi);
  for my $key (@keys) {
    my ($emin, $emax) = $stats{$key}{energy}->minmax;
    my $ymin = ($stats{$key}{mean} - $stats{$key}{prms})->min;
    my $ymax = ($stats{$key}{mean} + $stats{$key}{prms})->max;

    $elo = $emin unless defined $elo and $elo < $emin;
    $ehi = $emax unless defined $ehi and $ehi > $emax;
    $ylo = $ymin unless defined $ylo and $ylo < $ymin;
    $yhi = $ymax unless defined $yhi and $yhi > $ymax;
  }

  my $margin;

  $margin = ($ehi - $elo) / 10;
  $elo -= $margin;
  $ehi += $margin;

  $margin = ($yhi - $ylo) / 10;
  $ylo -= $margin;
  $yhi += $margin;

  pgsave();
#  pgsch($opts{charsize});
  pgenv($elo, $ehi, $ylo, $yhi, 0, 10);
  pglabel('energy (keV)', '', '');
  Chandra::Tools::Common::legend( text => \@legends,
				  points => \@symbols,
				  colors => \@colors,
				);
  for my $i (0..$#keys) {
    my $key = $keys[$i];
    my $symbol = $symbols[$i];
    my $color = $colors[$i];
    my $energy = $stats{$key}{energy}->float;

    # move the energy every so slightly so everything isn't overlapping
    $energy += ($ehi-$elo)*.01 * ($i - (@keys/2-.5));

    my $mean = $stats{$key}{mean}->float;
    my $prms = $stats{$key}{prms}->float;
    pgsave();
    pgsci($color);
    pgpt($energy->nelem, $energy->get_dataref, $mean->get_dataref, $symbol);
    pgerry($energy->nelem, $energy->get_dataref, ($mean-$prms)->get_dataref, ($mean+$prms)->get_dataref, 1);
    pgunsa();
  }
  pgunsa();

}
exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
