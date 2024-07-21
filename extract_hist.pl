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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> January 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use Lab;

use Getopt::Long;
my %default_opts = (
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'info!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

if (($opts{info} and @ARGV != 1) or (!$opts{info} and @ARGV != 5)) {
  die "usage: $0 (--info | v vsub u usub) binfile\n";
}

if ($opts{info}) {
  my $bin = shift;
  my $it = Lab::InBinFile->new($bin) or die;
  print "crsv\tvsub\tcrsu\tusub\ty1\ty2\tx1\tx2\n";
  print "N\tN\tN\tN\tN\tN\tN\tN\n";
  while (my ($crsv, $vsub, $crsu, $usub, $y1, $y2, $x1, $x2) =
	 $it->next_subtap(0)) {
    print join("\t",
	       $crsv, $vsub, $crsu, $usub, $y1, $y2, $x1, $x2,
	       ),"\n";
  }
}

# extact the histogram information for the requested subtap
else {
  my ($crsv, $vsub, $crsu, $usub, $bin) = @ARGV;
  my $it = Lab::InBinFile->new($bin) or die;
  while (
	 my ($tcrsv, $tvsub, $tcrsu, $tusub, $y1, $y2, $x1, $x2, $x, $y) = $it->next_subtap
	) {
    if (
	$tcrsv == $crsv and
	$tvsub == $vsub and
	$tcrsu == $crsu and
	$tusub == $usub
       ) {
      print "pha\tn\n";
      print "N\tN\n";
      $PDL::IO::Misc::colsep='';
      wcols "%d\t%d", $x, $y;
      exit;
    }
  }

  die "did not find requested subtap\n";
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
