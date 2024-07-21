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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> July 2005

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::GSL::RNG;
use Lab;
use Chandra::Tools::Common qw( read_bintbl_cols );

use Getopt::Long;
my %default_opts = (
		    nsamp => 1e4,
		    dev => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'nsamp=i', 'dev=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my $rng = PDL::GSL::RNG->new('default');

# read SAMP data, just a single chip's data for now
my ($mcp, $hrcfile) = (Lab::test_data('C-Ka'))[2,4];
$_ = $_->[0] for $mcp, $hrcfile;
print $hrcfile,"\n";
print $mcp,"\n";
my ($samp, $chip_id) = read_bintbl_cols("$Lab::EVTDIR/${hrcfile}_evt1_filt_spi.fits", qw( samp chip_id ), { extname => 'events' }) or die;

# filter on exposed plate
my $index = which($chip_id == Lab::mcp_to_chipid($mcp));
$_ = $_->index($index)->sever for $samp, $chip_id;

bin hist $samp;

my @n = (30, 100, 300, 1000);
my @trim = (5, 10, 15, 20);

dev $opts{dev}, 3, 2;

for my $n (@n) {
  my $sample = $samp->index(($rng->get_uniform($n,$opts{nsamp})*$n)->long)->qsort;
  my $median_b = $sample->medover;
  bin hist($median_b), { xtitle => sprintf("median (\\gm=%.3g, \\gs=%.3g), n=$n", ($median_b->stats)[0,1])};

  my $mean_b = $sample->average;
  bin hist($mean_b), { xtitle => sprintf("mean (\\gm=%.3g, \\gs=%.3g), n=$n", ($mean_b->stats)[0,1])};

  for my $t (@trim) {
    my $trimmed_mean_b = Lab::trim($sample, $t / 100, 1)->average;
    bin hist($trimmed_mean_b), { xtitle => sprintf("%d%% trimmed mean (\\gm=%.3g, \\gs=%.3g), n=%d", $t, ($trimmed_mean_b->stats)[0,1], $n)};
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
