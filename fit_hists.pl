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
use Lab;
use PDL;
use File::Temp;

use Getopt::Long;
my %default_opts = (
		    mincnts => 50,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'info!', 'mincnts=i',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my $bin = shift;

my $it = Lab::InBinFile->new($bin) or die;

my @cols = qw( ytap ysubtap xtap xsubtap y1 y2 x1 x2 sum );

if ($opts{info}) {

  print join("\t", @cols),"\n";
  print join("\t", ('N')x@cols),"\n";

  while (my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) =
	 $it->next_subtap) {
    print join("\t",
	       $ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $y->sum,
	       ),"\n";
  }
}


# do a single fit of the first dataset
else {

  push @cols, qw( g1fwhm g1pos g1ampl g2fwhm g2pos g2ampl );
  print join("\t", @cols),"\n";
  print join("\t", ('N')x@cols),"\n";

  while (
	 my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) = $it->next_subtap
	) {
    next unless $y->sum >= $opts{mincnts};

    my ($g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl) = fit_hist($x, $y);

    print join("\t",
	       $ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $y->sum,
	       $g1fwhm, $g1pos, $g1ampl, $g2fwhm, $g2pos, $g2ampl,
	      ),"\n";
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

  $cmd->print(<<EOP);
paramprompt off
()=load_dataset("$data")
source = ngauss1d[g1] + ngauss1d[g2]
g1.fwhm=30
g1.pos=@{[$mean/2]}
g1.ampl=@{[$y->sum / 4]}
g2.fwhm=20
g2.pos=@{[$mean]}
g1.ampl=@{[$y->sum / 2]}
ignore bins 1:1
ignore bins 256:256
show all
fit
EOP

  $cmd->close;

  my @lines = `sherpa --batch $cmd`;

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
