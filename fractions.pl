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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> June 2007

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
use Lab;

use Getopt::Long;
my %default_opts = (
		    rdb => 'hrcs_lab.rdb',
		    bindir => '/data/legs/rpete/data/hrcs_lab/analysis',
		    mincnts => 50,
		    dev => '/xs',
		    charsize => 2,
		    bgbsize => 4,
		    type => 'samp',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'rdb=s', 'datadir=s', 'mincnts=i', 'dev=s',
	   'charsize=f', 'maxmed!', 'maxmedn!', 'maxmedp!',
	   'dual!', 'type=s', 'bgbsize=i',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

# we'll toggle this on and off when we read the samp data
$opts{type} = 'pha' if $opts{maxmed} and $opts{dual};

if ($opts{maxmedp}) {
  warn "setting --bgbsize=1 for compatibility with --maxmedp\n";
  $opts{bgbsize} = 1;
  $opts{maxmed} = 1;
}

my @anodes = @ARGV ? @ARGV : @{ +( Lab::lines() )[0] };
my @energies = Lab::energies(@anodes);

dev $opts{dev}, ($opts{maxmed} ? (3, 3) : (3, 2));

# fractions in percent
my @fracs = (0.10, 0.50, 1.25, 96.25, 98.50, 99.60);

my @a = (@anodes, ($opts{maxmed} ? ('Background') : ()));

my %maxmedp; # if printing maxmed histograms

for my $i (0..$#a) {
  my $anode = $a[$i];

  # default values for the case of background plot
  my ($mcp, $hrc_file) = ([], ['merged_bg']);

  if ($anode ne 'Background') {
    (undef, undef, $mcp, undef, $hrc_file) = Lab::test_data($anode)
  }

  my ($n, $hists) = process($hrc_file, $mcp);
  my @n = @$n;
  my @hists = @$hists;

  if (!$opts{maxmed}) {
    for my $j (0..$#fracs) {
      my $frac = $fracs[$j];

      # channel at which this fraction is exceeded for each median
      my $n = $opts{type} eq 'samp' ? 512 : 256;
      my @lower = (0)x$n;

      for my $i (0..$#n) {
	next unless $n[$i];

=begin comment

      my $index =which($hists[$i]->cumusumover->double / $hists[$i]->sum <= $frac/100);
      $index->nelem or next;

      $lower[$i] = $index->at(-1);

=cut

	$lower[$i] = $hists[$i]->rld(sequence(long,$n))->pct($frac/100);
      }

      my $lower = float(\@lower);
      my $median = sequence(float, $n);

      my $i = which($lower > 0);
      points $median->index($i), $lower->index($i)/$median->index($i), { ytitle => sprintf("PHA\\d%.2f%%\\u / PHA\\dmedian\\u",$frac), xtitle => 'PHA\\dmedian\\u', ($j==1?(title=>$anode):()) };
    }
  }

  # plot distribution for maximum median
  if ($opts{maxmed}) {

    my $maxi = maximum_ind(float(\@n));
    my $x = sequence(long, 256);
    my $y = $hists[$maxi]->slice('0:255');

    if ($anode eq 'Background') {
      # rebin
#      my $n = $opts{type} eq 'samp' ? 512 : 256;
      my $n = 256;
      my $vals = $y->rld(sequence(long, $n));
      ($x, $y) = hist $vals->where($vals<256), -0.5, 255.5, $opts{bgbsize};
    }

    # get SAMP data if we're doing a dual plot
    my ($x2, $y2);
    if ($opts{dual}) {
      $opts{type} = 'samp';
      my ($n, $hists) = process($hrc_file, $mcp);
      my @n = @$n;
      my @hists = @$hists;
      my $maxi = maximum_ind(float(\@n));
      my $x = sequence(long, 256);
      my $y = $hists[$maxi]->slice('0:255');

      if ($anode eq 'Background') {
	# rebin
	#      my $n = $opts{type} eq 'samp' ? 512 : 256;
	my $n = 256;
	my $vals = $y->rld(sequence(long, $n));
	($x, $y) = hist $vals->where($vals<256), -0.5, 255.5, $opts{bgbsize};
      }
      $opts{type} = 'pha';
      ($x2, $y2) = ($x, $y);
    }

    my $ymax = $y->max * 1.1;
    $ymax = $y2->max*1.1 if $opts{dual} and $y2->max > $y->max;

    my $hardlw = 2;
    my $hardch = 2.2;

    pgsave();

    pgslw($hardlw);
    pgsch($hardch);

    my $xtitle = '';

    # x axis label on bottom row of 3x3 grid
    if ($i > 5) {
      $xtitle = uc($opts{type});
      $xtitle = 'channel' if $opts{dual};
    }

    my $xopts = 'BCT';
    my $yopts = $xopts;

    $xopts .= 'N' if $i > 5; # bottom row gets numeric labels
    $yopts .= 'N' if $opts{maxmedn}; # don't put numbers on y axis by default

    pgenv(0, 255, 0, $ymax, 0, -2);
    pgbox($xopts, 0, 0, $yopts, 0, 0);
    pglabel($xtitle, '', $anode);

    pgsci(4);
    pgbin($x->nelem, $x->float->get_dataref, $y->float->get_dataref, 1);

    if ($opts{maxmedp}) {
      $maxmedp{$opts{type}} = $x;
      $maxmedp{$anode} = $y;
    }

#    bin ($x, $y, { title => $anode, xtitle => $xtitle, yrange=>[0,$ymax], color => 'blue', hardlw => $hardlw, hardch => $hardch, }, );
    if ($opts{dual}) {
      pgsci(2);
      pgbin($x2->nelem, $x2->float->get_dataref, $y2->float->get_dataref, 1);
#      hold;
#      bin $x2, $y2, { color => 'red', hardlw => $hardlw, hardch => $hardch, };
#      release;
    }
#    line($hists[$maxi], { title => $anode });
#    if ($anode eq 'Al-Ka') {
#      line($hists[90+$_], { title => "$anode, median = ".(90+$_) }) for 0..11;
#    }

    pgunsa();
  }

  1;
}

if ($opts{maxmedp}) {
  print join("\t", $opts{type}, @a), "\n";
  print join("\t", ('N')x(@a+1)),"\n";
  my $fmt = join "\t", ('%d')x(@a+1);
  $PDL::IO::Misc::colsep = '';
  wcols $fmt, @maxmedp{$opts{type}, @a}, *STDOUT{IO};
}

exit 0;

sub process {
  my ($files, $mcp) = @_;

  $Lab::TYPE = $opts{type};

  my $n = $opts{type} eq 'samp' ? 512 : 256;

  my @n = (0)x$n;
  my @hists = map { zeroes(long, $n) } 1..$n;

  my $mcp_filter = @$mcp;

  for (0..$#{$files}) {

    my $file = $files->[$_];
    my $mcp = $mcp->[$_] if $mcp_filter;

    my $bin = "$opts{bindir}/${file}_$opts{type}.bin";

    print STDERR "Processing $bin...";
    STDERR->flush;

    my ($rawy_min, $rawy_max) = Lab::rawy_limits_mcp($mcp) if $mcp_filter;

    my $it = Lab::InBinFile->new($bin) or die;

    while (my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) = $it->next_subtap) {

      next unless !$mcp_filter or ($y1 >= $rawy_min && $y2 <= $rawy_max);

      my $sum = $y->sum;
      $sum >= $opts{mincnts} or next;

      my $vals = $y->rld($x);
      my $median = rint($vals->median)->at;

      $n[$median]++;
      $hists[$median] += $y;
    }

    print STDERR "done\n";
  }

  return \(@n, @hists);
}

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
