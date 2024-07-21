#! /usr/bin/perl -w
use strict;

=head1 NAME

height_vs_energy.pl - plots of PHA, SAMP and SPI versus log10 energy

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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> June 2008

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use Chandra::Tools::Common qw( read_bintbl_cols );
use Lab;
use IO::File;

use Getopt::Long;
my %default_opts = (
		    analdir => $Lab::ANALDIR,
		    dev => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'analdir=s', 'dev=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my @anodes = @ARGV;

if (!@anodes) {
  my ($lines) = Lab::lines();
  @anodes = @$lines;
}

my @energies = Lab::energies(@anodes);
my $e = pdl(\@energies)->log10;

dev $opts{dev}, 3, scalar(@anodes);

my %mean = ( samp => [], spi => [], pha => [] );
my %rms = ( samp => [], spi => [], pha => [] );

# for each anode, read fits file pha/samp/spi, filter on mcp, plot them
for my $l (@anodes) {

  my ($mcp, $hrcfile) = (Lab::test_data($l))[2,4];

  my $pha = long [];
  my $samp = double [];
  my $spi = $samp->copy;

  for my $j (0..$#{$mcp}) {

    my ($mcp, $hrcfile) = ($mcp->[$j], $hrcfile->[$j]);

    my $chip = Lab::mcp_to_chipid($mcp);

    my $fits = "$opts{analdir}/${hrcfile}_evt1_filt_spi.fits";

    print STDERR "reading $fits...";
    STDERR->flush;

     my ($p, $sa, $sp) = read_bintbl_cols($fits, qw( pha samp spi ), { extname => 'events', rfilter => "chip_id == $chip"});
     print STDERR " done\n";
     $pha = append($pha, $p);
     $samp = append($samp, $sa);
     $spi = append($spi, $sp);
  }

  $_ = Lab::trim($_, .05) for $pha, $samp, $spi;

  bin hist( $pha), { xtitle => 'PHA', title => $l };
  bin hist( $samp), { xtitle => 'SAMP', title => $l };
  bin hist( $spi), { xtitle => 'SPI', title => $l };

  # calculate mean and std dev, store them
  my ($m, $r);

  ($m, $r) = $pha->double->stats;
  push @{$mean{pha}}, $m->at;
  push @{$rms{pha}}, $r->at;

  ($m, $r) = $samp->double->stats;
  push @{$mean{samp}}, $m->at;
  push @{$rms{samp}}, $r->at;

  ($m, $r) = $spi->double->where($spi==$spi)->stats;
  push @{$mean{spi}}, $m->at;
  push @{$rms{spi}}, $r->at;
}


# plot mean of pha/samp/spi vs energy with error bars
points $e, pdl($mean{pha}), { ytitle => 'average PHA', xtitle => 'log E', border => 1 };
points $e, pdl($mean{samp}), { ytitle => 'average SAMP', xtitle => 'log E', border => 1 };
points $e, pdl($mean{spi}), { ytitle => 'average SPI', xtitle => 'log E', border => 1 };

points $e, pdl($rms{pha}), { ytitle => 'rms PHA', xtitle => 'log E', border => 1 };
points $e, pdl($rms{samp}), { ytitle => 'rms SAMP', xtitle => 'log E', border => 1 };
points $e, pdl($rms{spi}), { ytitle => 'rms SPI', xtitle => 'log E', border => 1 };

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
