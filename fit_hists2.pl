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
use PDL::Graphics::PGPLOT::Window;
use PGPLOT;
use MyRDB;
use Data::Dumper;
use File::Temp;
use Lab;
use PDL::Fit::Polynomial;

use Getopt::Long;
my %default_opts = (
		    rdb => 'hrcs_lab.rdb',
		    bindir => '/data/legs/rpete/data/hrcs_lab/analysis',
		    mincnts => 50,
		    minsubtaps => 10,
		    dev => '/xs',
		    charsize => 2,
		    fitfits => 1,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'rdb=s', 'datadir=s', 'mincnts=i', 'dev=s', 'minsubtaps=i',
	   'revcolors!', 'colors!', 'charsize=f', 'fitfits!', 'samp!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my @anodes = @ARGV ? @ARGV : @{ +( Lab::lines() )[0] };
my @energies = Lab::energies(@anodes);

my $nz = 256;
my $name = 'PHA';
if ($opts{samp}) {
  $nz = 512;
  $name = 'SAMP';
  $Lab::SAMP = 1;
}

my $dev = PDL::Graphics::PGPLOT::Window->new(Device => $opts{dev},
					     NXPanel => 4, NYPanel => 3,
					     ) or die;

my @cols = qw( median g1fwhm g1pos g1ampl g2fwhm g2pos g2ampl );

print join("\t", @cols),"\n";
print join("\t", ('N')x@cols),"\n";

my (@u2_median, @u2_median_err,
    @s2_median, @s2_median_err,
    @u1_u2, @u1_u2_err,
    @n1_n2, @n1_n2_err,
   );

for my $anode (@anodes) {

  my (undef, undef, $MCP, undef, $HRC_file) = Lab::test_data($anode);

  my @n = (0)x$nz;
  my @hists = map { zeroes(long, $nz) } 1..$nz;

  for (0..$#{$HRC_file}) {

    my $hrc_file = $HRC_file->[$_];
    my $mcp = $MCP->[$_];

    my $bin = $opts{bindir} . '/' . $hrc_file . ($opts{samp} ? '_samp' : '') . '.bin';

    print STDERR "Processing $bin...";
    STDERR->flush;

    my ($rawy_min, $rawy_max) = Lab::rawy_limits_mcp($mcp);

    my $it = Lab::InBinFile->new($bin) or die;

    while (my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) = $it->next_subtap) {

      next if $y1 < $rawy_min or $y2 > $rawy_max;

      my $sum = $y->sum;
      $sum >= $opts{mincnts} or next;

      my $vals = $y->rld($x); # PHA or SAMP values for events in this subtap
      my $median = rint($vals->median)->at;

      $n[$median]++;
      $hists[$median] += $y;
    }

    print STDERR "done\n";
  }

  my (@nsubtaps, @ncounts, @median, @g1fwhm, @g2fwhm, @g1pos, @g2pos, @g1ampl, @g2ampl);

  for my $i (0..$#n) {
    $n[$i] >= $opts{minsubtaps} or next;

    my ($x, $y) = (sequence(long, $nz), $hists[$i]);

#    line $x, $y, { title => "$anode, median = $i" };

    my ($g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl) = fit_hist($x,$y);

    print join("\t",
	       $i, $g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl),"\n";

    push @nsubtaps, $n[$i];
    push @ncounts, $y->sum;
    push @median, $i;
    push @g1fwhm, $g1fwhm;
    push @g2fwhm, $g2fwhm;
    push @g1pos, $g1pos;
    push @g2pos, $g2pos;
    push @g1ampl, $g1ampl;
    push @g2ampl, $g2ampl;

  }

  my $nsubtaps = pdl \@nsubtaps;
  my $ncounts = pdl \@ncounts;
  my $median = pdl \@median;
  my $g1fwhm = pdl \@g1fwhm;
  my $g2fwhm = pdl \@g2fwhm;
  my $g1pos = pdl \@g1pos;
  my $g2pos = pdl \@g2pos;
  my $g1ampl = pdl \@g1ampl;
  my $g2ampl = pdl \@g2ampl;

  $dev->points($median, $g1pos, { title => $anode, xtitle => "median $name", ytitle => '\\gm\\d1\\u' });
  $dev->points( $median, $g2pos, { title => $anode, xtitle => "median $name", ytitle => '\\gm\\d2\\u' });
  $dev->points( $median, $g1pos/$median, { title => $anode, xtitle => "median $name", ytitle => "\\gm\\d1\\u / median $name" });
  $dev->points( $median, $g2pos/$median, { title => $anode, xtitle => "median $name", ytitle => "\\gm\\d2\\u / median $name" });
  $dev->points( $median, $g1ampl / $ncounts, { title => $anode, xtitle => "median $name", ytitle => 'N\\d1\\u' });
  $dev->points( $median, $g2ampl / $ncounts, { title => $anode, xtitle => "median $name", ytitle => 'N\\d2\\u' });
  $dev->points( $median, $g1fwhm/2.354, { title => $anode, xtitle => "median $name", ytitle => '\\gs\\d1\\u' });
  $dev->points( $median, $g2fwhm/2.354, { title => $anode, xtitle => "median $name", ytitle => '\\gs\\d2\\u' });
  $dev->points( $median, $g2fwhm/2.354 / $median, { title => $anode, xtitle => "median $name", ytitle => "\\gs\\d2\\u / median $name" });
  $dev->points( $median, $g1pos/$g2pos, { title => $anode, xtitle => "median $name", ytitle => '\\gm\\d1\\u / \\gm\\d2\\u' });
  $dev->points( $median, $g1ampl/$g2ampl, { title => $anode, xtitle => "median $name", ytitle => 'N\\d1\\u / N\\d2\\u' });
  $dev->points( $median, $g1fwhm/$g2fwhm, { title => $anode, xtitle => "median $name", ytitle => '\\gs\\d1\\u / \\gs\\d2\\u' });

  my $u2_median = sum($g2pos / $median * $ncounts / $ncounts->sum);
  my $u2_median_err = +(stats($g2pos / $median))[1]->at;

  my $s2_median = sum($g2fwhm / 2.354 / $median * $ncounts / $ncounts->sum);
  my $s2_median_err = +(stats($g2fwhm / 2.354 / $median))[1]->at;

  my $u1_u2 = sum($g1pos / $g2pos * $ncounts / $ncounts->sum);
  my $u1_u2_err = +(stats($g1pos / $g2pos))[1]->at;

  my $n1_n2 = sum($g1ampl / $g2ampl * $ncounts / $ncounts->sum);
  my $n1_n2_err = +(stats($g1ampl / $g2ampl))[1]->at;

  push @u2_median, $u2_median;
  push @u2_median_err, $u2_median_err;
  push @s2_median, $s2_median;
  push @s2_median_err, $s2_median_err;
  push @u1_u2, $u1_u2;
  push @u1_u2_err, $u1_u2_err;
  push @n1_n2, $n1_n2;
  push @n1_n2_err, $n1_n2_err;
}

#pgsubp(3,1);

if ($opts{fitfits}) {
   my $energy = pdl(\@energies) * 1e-3; # convert to keV
   my $u2_median = pdl \@u2_median;
   my $u2_median_err = pdl \@u2_median_err;
   my $s2_median = pdl \@s2_median;
   my $s2_median_err = pdl \@s2_median_err;
   my $u1_u2 = pdl \@u1_u2;
   my $u1_u2_err = pdl \@u1_u2_err;
   my $n1_n2 = pdl \@n1_n2;
   my $n1_n2_err = pdl \@n1_n2_err;

   my $x = $energy->log10;

   my ($xlo, $xhi) = $x->minmax;
   my $xdiff = 0.1 * ($xhi-$xlo);
   $xlo -= $xdiff; $xhi += $xdiff;

   for (
        [$u2_median, $u2_median_err, "\\gm\\d2\\u / median $name"],
        [$s2_median, $s2_median_err, "\\gs\\d2\\u / median $name"],
        [$u1_u2, $u1_u2_err, '\\gm\\d1\\u / \\gm\\d2\\u'],
        [$n1_n2, $n1_n2_err, 'N\\d1\\u / N\\d2\\u'],
       ) {

     my ($y, $yerr) = ($_->[0], $_->[1]);
     my $ytitle = $_->[2];

     my $yhi = +($y + $yerr)->max;
     my $ylo = +($y - $yerr)->min;

     my $ydiff = 0.1 * ($yhi-$ylo);
     $ylo -= $ydiff; $yhi += $ydiff;

     $dev->points($x, $y, { axis => 'logx', xrange => [$xlo,$xhi], yrange => [$ylo,$yhi], xtitle => 'energy (keV)', ytitle => $ytitle });
     $dev->hold;

     $dev->errb($x, $y, $yerr);
     $dev->hold;

     my ($yfit, $coeffs) = fitpoly1d($x, $y, 2, { Weights => 1/$yerr/$yerr });
     $dev->line(pdl($xlo,$xhi), pdl($xlo, $xhi)*$coeffs->at(1)+$coeffs->at(0));

     $dev->release;
   }

}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub fit_hist {
  my ($x, $y) = @_;

  my $data = File::Temp->new
    or die "could not create temporary data file: $!";;

  wcols($x, $y, $data);
  $data->close;

  my $cmd = File::Temp->new
    or die "could not create temporary command file: $!";

  # use mean as main peak initial value.
  my $mean = sum($x * $y) / $y->sum;

  my $script = <<EOP;
paramprompt off
()=load_dataset("$data")
source = ngauss1d[g1] + ngauss1d[g2]
g1.fwhm=50
g1.pos=@{[$mean/2]}
g1.ampl=@{[0.2 * $y->sum]}
g1.ampl.max=@{[0.6 * $y->sum]}
g2.fwhm=20
g2.pos=@{[$mean]}
g2.ampl=@{[0.8 * $y->sum]}
ignore bins 1:1
ignore bins $nz:$nz
show all
fit
EOP

#  print $script;

  $cmd->print($script);

  $cmd->close;

  my @lines = `sherpa --batch $cmd`;
  $? and die;

  my $real = qr/-?\d+\.?\d*/;
  my ($g1fwhm) = grep { /g1\.fwhm\s+$real\s+/ } @lines;
  my ($g1pos) = grep { /g1\.pos\s+$real\s+/ } @lines;
  my ($g1ampl) = grep { /g1\.ampl\s+$real\s+/ } @lines;
  my ($g2fwhm) = grep { /g2\.fwhm\s+$real\s+/ } @lines;
  my ($g2pos) = grep { /g2\.pos\s+$real\s+/ } @lines;
  my ($g2ampl) = grep { /g2\.ampl\s+$real\s+/ } @lines;

#  my ($g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl) = @lines[$#lines-6..$#lines-1];

  ($_) = /\s+($real)\s+/ for $g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl;
  return $g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl;
}

sub _version {
  print $version,"\n";
  exit 0;
}
