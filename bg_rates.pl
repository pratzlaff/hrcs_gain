#! /usr/bin/perl -w
use strict;

=head1 NAME

bg_extract.pl - whee

=head1 SYNOPSIS

cp template newprog

=head1 DESCRIPTION

Extract background data from unexposed plates in the HRC-S flat field
lab data.

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> February 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use MyRDB;
use PDL;
use PDL::Graphics::PGPLOT;
use Chandra::Tools::Common qw( read_bintbl_cols );
use IO::Handle;
use Data::Dumper;
use Lab;

use Getopt::Long;
my %default_opts = (
		    evtdir => $Lab::EVTDIR,
		    bgdir => $Lab::BGDIR,
		    config => $Lab::TESTFILE,
		    filt => 1,
		    binsize => 10,
		    dev => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'evtdir=s', 'bgdir=s', 'config=s', 'filt!', 'binsize=f', 'dev=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV and die "usage: $0 [options]\n";

$Lab::TESTFILE = $opts{config};

my (undef, undef, $mcp, $evt_time, $hrc_file, $bg_time, $bg_file)
   = Lab::test_data() or die "problem reading $opts{config}";

# pdls of X,Y event positions, exposure times for each chip_id
my (@rawx, @rawy, @exptime);

# initialze our data
for (0..2) {
  $rawx[$_] = long [];
  $rawy[$_] = $rawx[$_]->copy;
  $exptime[$_] = 0;
}

my ($rawx, $rawy, $chip_id);

my $suffix = '_evt1.fits';
$suffix = '_evt1_filt.fits' if $opts{filt};

for my $i (0..$#{$hrc_file}) {

  my $mcp = $mcp->[$i];
  my $hrc_file = $hrc_file->[$i];
  my $bg_file = $bg_file->[$i];
  my $bg_time = $bg_time->[$i];
  my $evt_time = $evt_time->[$i];

  my $evtfile = $opts{evtdir}. '/' . $hrc_file . $suffix;
  my $bgfile = $opts{bgdir}. '/' . $bg_file . $suffix;

  # read bg data
  print STDERR "reading $bgfile...";
  STDERR->flush;
  ($rawx, $rawy, $chip_id) = read_bintbl_cols($bgfile, 'rawx', 'rawy', 'chip_id', { status => 1, extname => 'events' }) or die;
  print STDERR " done\n";

  for (0..2) {
    my $ind = which($chip_id == $_+1);
    $rawx[$_] = $rawx[$_]->append($rawx->index($ind));
    $rawy[$_] = $rawy[$_]->append($rawy->index($ind));
    $exptime[$_] += $bg_time;
  }

  # do same for evt file, skipping current active MCP
  print STDERR "reading $evtfile...";
  STDERR->flush;
  ($rawx, $rawy, $chip_id) = read_bintbl_cols($evtfile, 'rawx', 'rawy', 'chip_id', { status => 1, extname => 'events' }) or die;
  print STDERR " done\n";

  my $chip_id_active = Lab::mcp_to_chipid($mcp);;
  for (0..2) {
    next if $_+1 == $chip_id_active;
    my $ind = which($chip_id == $_+1);
    $rawx[$_] = $rawx[$_]->append($rawx->index($ind));
    $rawy[$_] = $rawy[$_]->append($rawy->index($ind));
    $exptime[$_] += $evt_time;
  }

}

dev $opts{dev}, 2, 3;

for (0..2) {
  my $rawx = $rawx[$_];
  my $rawy = $rawx[$_];
  my $exptime = $exptime[$_];
  my $chip_id = $_ + 1;

  my ($xmin, $xmax) = $rawx->minmax;
  my ($ymin, $ymax) = $rawy->minmax;

  my ($histx, $histy);

  # rawx histogram, counts / sec / pixel
  ($histx, $histy) = hist($rawx, $xmin, $xmax, $opts{binsize});
  $histy = float($histy) / $opts{binsize} / ($ymax-$ymin) / $opts{binsize} / $exptime;
  bin $histx, $histy, { xtitle => 'rawx', ytitle => 'cnts / sec / pixel\\u2\\d', title => "chip_id $chip_id" };

  # rawy histogram, counts / sec / pixel
  ($histx, $histy) = hist($rawy, $ymin, $ymax, $opts{binsize});
  $histy = float($histy) / $opts{binsize} / ($xmax-$xmin) / $opts{binsize} / $exptime;
  bin $histx, $histy, { xtitle => 'rawy', ytitle => 'cnts / sec / pixel\\u2\\d', title => "chip_id $chip_id" };

}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
