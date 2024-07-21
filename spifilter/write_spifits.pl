#! /usr/bin/perl -w
use strict;

=head1 NAME

write_spifits.pl - create FITS file used by addspi, addpilab

=head1 SYNOPSIS

perl write_spifits.pl spimeanfits2.rdb tgain1.oldstyle outfile.fits

=head1 DESCRIPTION

Reads Brad's parameter and tgain files, generates a FITS file suitable
to be used by addspi

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> June 2018

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Carp;
use PDL;
use Astro::FITS::CFITSIO qw( TDOUBLE BINARY_TBL );

use Getopt::Long;
my %default_opts = (
		    );

my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}

$opts{help} and _help();
$opts{version} and _version();

@ARGV == 3 or die "Usage: $0 [options] spifits tgain outfile\ntry --help for more information\n";

my ($spifits, $tgain, $outfile) = @ARGV;

my ($S, $b, $m) = read_spirdb($spifits);
my ($dates, $factors) = read_tgain($tgain);
write_spifits($outfile, $S, $b, $m, $dates, $factors);

exit 0;

sub _help {
  exec('perldoc', '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}


sub check_status {
  my $s = shift;
  if ($s != 0) {
    my $txt;
    Astro::FITS::CFITSIO::fits_get_errstatus($s,$txt);
    carp "CFITSIO error: $txt";
    return 0;
  }

  return 1;
}

# create the output file
sub write_spifits {

  # write everything
  #wfits cat($S, $b, $m)->mv(2,0)->float, $file;

  my ($file, $S, $b, $m, $dates, $factors) = @_;

  print STDERR "writing $file...";
  STDERR->flush;

  my $status=0;
  my $fptr = Astro::FITS::CFITSIO::create_file('!'.$file, $status);
  check_status($status) or die "error creating $file\n";

  # write "stacked" S_i, b_i and m_i as primary hdu
  my $img = cat($S, $b, $m)->mv(2,0)->double;
  $fptr->create_img(-64, $img->ndims, [$img->dims], $status);
  $fptr->write_pix(TDOUBLE, [(1)x$img->ndims], $img->nelem, $img->get_dataref, $status);
  $fptr->write_chksum($status);
  check_status($status) or die "error creating primary HDU in $file\n";

  # add binary table for gaindrop factors
  $fptr->create_tbl(BINARY_TBL, 0, 2,
		    ['date', 'dropcorr'],
		    ['1D', $factors->getdim(1).'D'],
		    undef,
		    'gaindrop',
		    $status);
  check_status($status) or die "error creating gaindrop HDU in $file\n";
  $_ = $_->double for $dates, $factors;
  $fptr->write_col(TDOUBLE, 1, 1, 1, $dates->nelem, $dates->get_dataref, $status);
  $fptr->write_col(TDOUBLE, 2, 1, 1, $factors->nelem, $factors->mv(0,1)->get_dataref, $status);
  $fptr->write_chksum($status);

  $fptr->close_file($status);
  check_status($status) or die "error closing $file\n";

  print STDERR "done\n";

}

sub read_spirdb {
  my $rdbfile = shift;

  my $S = ones(48, 576); # FIXME: are these dimensions always going to work?
  my $b = zeros($S);
  my $m = $S->copy;

  my ($v, $vsub, $u, $usub, $S_i, $b_i, $m_i) = rcols($rdbfile, { lines => '2:' });

  my %out = ( S => $S, b => $b, m => $m );
  my %in = ( S => $S_i, b => $b_i, m => $m_i );

  # FIXME: figure out what's going on with the indexing here
  my $yi = ($v-0)*3 + $vsub;
  my $xi = ($u-0)*3 + $usub;

  my $tmp;

  ($tmp = $out{$_}->index2d($xi, $yi)->flat) .= $in{$_} for keys %out;

  return @out{ qw/ S b m / };
}

sub read_tgain {
  open(my $fh, '<', $_[0]) or die $!;
  <$fh>; # discard first line
  my @times = split ' ', <$fh>;
  my $timegrid = pdl(@times[2..$#times]);
  my ($rawy, @tgain) = rcols $fh;
  my $tgain = cat @tgain;
  return 1998 + ($timegrid-50814)/365.2524, $tgain->xchg(0,1);
}
