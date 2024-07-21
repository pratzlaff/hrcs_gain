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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> May 2008

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use Lab;
use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Fit::Polynomial;
use FindBin;
use Config;
use Carp;
use MyRDB;

use Getopt::Long;
my %default_opts = (
		    dev => '/xs',
		    n => 500,
		    rdbdir => $Lab::ANALDIR,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'dev=s', 'n=i', 'rdbdir=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my @anodes = ();

# plot individual histograms for subtaps with IQR and RMS in this range
my ($il, $ih, $rl, $rh);

for my $n (scalar @ARGV) {
  $n == 0 and last;
  $n == 1 and @anodes = @ARGV, last;
  $n == 5 and @anodes = shift, ($il, $ih, $rl, $rh)=@ARGV, last;
}

my ($line, $energy, $HRC_file) = (Lab::test_data(@anodes))[0,1,4];

dev $opts{dev}, 6, 3;#, { aspectratio => 1/2 };

my (@intercept, @slope);

for my $i (0..$#{$line}) {
  my ($line, $energy, $hrc_file) = ($line->[$i], $energy->[$i], $HRC_file->[$i]);

  my $file = "$opts{rdbdir}/$hrc_file.rdb";
  -f $file or die $file;
  my ($v, $vsub, $u, $usub, $n, $srms, $siqr) = MyRDB::rdb_cols($file, qw( crsv vsub crsu usub n srms siqr ));

  $_ = pdl $_ for $v, $vsub, $u, $usub, $n, $srms, $siqr;
  my $i = which($n > $opts{n});
  $i->nelem or die "$opts{n} : $file";

  $_ = $_->index($i) for $v, $vsub, $u, $usub, $n, $srms, $siqr;

  points $srms, $siqr, {
			border => 1,
			title => $line . ' : ' . $hrc_file,
			xtitle => 'rms',
			ytitle => 'IQR',
		       };
  hold;
  my ($yfit, $coeff) = fitpoly1d($srms, $siqr, 2);
  line $srms, $yfit;
  release;

  push @intercept, $coeff->at(0);
  push @slope, $coeff->at(1);

  if (defined $il) {
    my $iqr = $siqr;
    my $rms = $srms;
    my $i = which(
		  ($iqr >= $il) & ($iqr <= $ih) &
		  ($rms >= $rl) & ($rms <= $rh)
		 );
    printf STDERR "%d subtaps found with IQR in [$il, $ih] and RMS in [$rl, $rh]\n", $i->nelem;
    local $PDL::IO::Misc::colsep = '';
    print join("\t", qw( crsv vsub crsu usub )),"\n";
    wcols "%d\t%d\t%d\t%d", $v->index($i), $vsub->index($i), $u->index($i), $usub->index($i);
  }
}

my $e = pdl($energy)->log10;
my $intercept = pdl \@intercept;
my $slope = pdl \@slope;
points $e, $intercept, {
			border => 1,
			xtitle => 'log \fiE',
			ytitle => 'intercept',
		       };
points $e, $slope, {
		    border => 1,
		    xtitle => 'log \fiE',
		    ytitle => 'slope',
		   };

printf "mean intercept = %f\n", $intercept->avg;
printf "mean slope = %f\n", $slope->avg;
printf "10%% trimmed mean intercept = %f\n", Lab::trim($intercept,.1)->avg;
printf "10%% trimmed mean slope = %f\n", Lab::trim($slope,.1)->avg;
printf "20%% trimmed mean intercept = %f\n", Lab::trim($intercept,.2)->avg;
printf "20%% trimmed mean slope = %f\n", Lab::trim($slope,.2)->avg;



exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
