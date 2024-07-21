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
use Astro::FITS::CFITSIO;
use Chandra::Tools::Common qw( check_status );
use IO::Handle;
use Lab;

use Getopt::Long;
my %default_opts = (
		    evt1dir => '/data/legs/rpete/data/hrcs_lab/evt1',
		    config => 'hrcs_lab.rdb',
		    omit => 1,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'config=s', 'omit!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV ==1 or die "usage: $0 [options] outfile\n";

my $args = "@ARGV"; # preserve original argument list for HISTORY entry

my $outfile = shift;

my ($MCP, $HRC_file) = MyRDB::rdb_cols($opts{config}, qw( MCP HRC_file ))
  or die "problem reading $opts{config}";

my @omit = Lab::merged_bg_omit();
my %omit; @omit{@omit} = ();

# fits_copy_file to start the output file, move to the events extension,
# fits_select_rows to filter on chip_id, then open rem

my $outfptr;
my $outstatus = 0;

for my $i (0..$#{$HRC_file}) {

  if ($opts{omit} and exists $omit{$HRC_file->[$i]}) {
    print STDERR "omitting $HRC_file->[$i]\n";
    next;
  }

  my $infile = $opts{evt1dir}. '/' . $HRC_file->[$i] . '_evt1.fits';

  print STDERR "processing $infile...";
  STDERR->flush;

  -f $infile or die "$infile not found";

  my $instatus = 0;
  my $infptr = Astro::FITS::CFITSIO::open_file($infile,
					       Astro::FITS::CFITSIO::READONLY(),
					       $instatus);

  # create output file if not done already
  if (!$outfptr) {
    $outfptr = Astro::FITS::CFITSIO::create_file('!'.$outfile, $outstatus);
    check_status($outstatus) or die "could not create $outfile";

    # copy all of first input file to output file
    $infptr->copy_file($outfptr, 1, 1, 1, $outstatus);
    check_status($outstatus) or die "error copying first input file";

    # move to the events extension, delete its data
    $outfptr->movnam_hdu(Astro::FITS::CFITSIO::BINARY_TBL(), 'events', 0, $outstatus);
    check_status($outstatus) or die "could not find events HDU in $infile";

    my $nrows;
    $outfptr->get_num_rows($nrows, $outstatus);
    $outfptr->delete_rows(1, $nrows, $outstatus);
    check_status($outstatus) or die "could not delete rows from output file";

    # add another history entry
    $outfptr->write_history("$0 $args", $outstatus);
  }

  # move to the input events extension
  $infptr->movnam_hdu(Astro::FITS::CFITSIO::BINARY_TBL(), 'events', 0, $instatus);
  check_status($instatus) or die "could not find events HDU in $infile";

  # copy select rows to output file
  my $expr;
  for my $mcp ($MCP->[$i]+0) {
    $mcp == 1 and $expr='chip_id != 1', last;
    $mcp == 0 and $expr='chip_id != 2', last;
    $mcp == -1 and $expr='chip_id != 3', last;
    die "unrecognized MCP = $mcp";
  }

  $infptr->select_rows($outfptr, $expr, $outstatus);
  check_status($outstatus) or die "error selecting rows from $infile";

  $infptr->close_file($instatus);
  check_status($instatus) or warn("nonzero status for $infile");

  print STDERR "done\n";
  STDERR->flush;
}

# update the DATASUM and CHECKSUM header entries
$outfptr->write_chksum($outstatus);

$outfptr->close_file($outstatus);
check_status($outstatus) or warn("nonzero status for $outfile");

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
