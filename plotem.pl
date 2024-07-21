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

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> February 2007

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use PDL;
use PDL::Graphics::PGPLOT;
use MyRDB;
use Data::Dumper;

use Getopt::Long;
my %default_opts = (
		    rdb => 'hrcs_lab.rdb',
		    datadir => '/data/legs/rpete/data/hrcs_lab/analysis',
		    mincnts => 400,
		    dev => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'rdb=s', 'datadir=s', 'mincnts=i', 'dev=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

my ($line, $energy, $hrc_file) = MyRDB::rdb_cols($opts{rdb}, qw( line energy HRC_file )) or die;

my @anodes = @ARGV ? @ARGV : @$line;

# remove duplicates
my %anodes;
@anodes{@anodes} = ();
@anodes = keys %anodes;

# ensure each exists in our rdb file
my %lines;
@lines{@$line} = @$energy;

for (@anodes) {
  exists $lines{$_} or die "no entries in $opts{rdb} for anode '$_'";
}

# order by energy
@anodes = sort { $lines{$a} <=> $lines{$b} } @anodes;

dev $opts{dev}, 2, 1;

for my $anode (@anodes) {

  my @i = grep { $line->[$_] eq $anode } 0..$#{$line};

  # energy (eV) for this anode
  my $energy = $energy->[$i[0]];

  my @hrc_file = @{$hrc_file}[@i];

  my $n = long [];
  my $mean = float [];
  my $median = $mean->copy;
  my $gsigma = $mean->copy;
  my $gmean = $gsigma->copy;
  my $g1pos = $gmean->copy;
  my $g2pos = $g1pos->copy;

  for my $hrc_file (@hrc_file) {
    my $rdb = $opts{datadir} . '/' . $hrc_file . '.rdb';
    my $rdb_gfits = $opts{datadir} . '/' . $hrc_file . '_gfits.rdb';
    for ($rdb, $rdb_gfits) {
      -f $_ or die "'$_' does not exist";
    }

    my ($tn, $tmean, $tmedian, $tgsigma, $tgmean)
      = MyRDB::rdb_cols( $rdb, qw( n mean median gsigma gmean ) )
	or die;

    my ($tsum, $tg1pos, $tg2pos)
      = MyRDB::rdb_cols( $rdb_gfits, qw( sum g1pos g2pos ) )
	or die;

    $_ = long $_ for $tn, $tsum;
    $_ = float $_ for $tmean, $tmedian, $tgsigma, $tgmean, $tg1pos, $tg2pos;

    my $i1 = which($tn >= $opts{mincnts});
    my $i2 = which($tsum >= $opts{mincnts});

    next unless $i1->nelem;

    $i1->nelem == $i2->nelem or die $hrc_file;

    # ensure we're looking at the same subtaps in both files
    my $sum_diff = $tn->index($i1) - $tsum->index($i2);
    $sum_diff->abs->max == 0 or die $hrc_file;

    $n = $n->append($tn->index($i1));
    $mean = $mean->append($tmean->index($i1));
    $median = $median->append($tmedian->index($i1));
    $gsigma = $gsigma->append($tgsigma->index($i1));
    $gmean = $gmean->append($tgmean->index($i1));

    $g1pos = $g1pos->append($tg1pos->index($i2));
    $g2pos = $g2pos->append($tg2pos->index($i2));
  }

  points $mean, $median, { xtitle => 'mean PHA', ytitle => 'median PHA', title => $anode, xrange => [20,200], yrange => [20,200], symbol => -1 };
  points $g1pos, $g2pos, { xtitle => '\\gm\\d1\\u', ytitle => '\\gm\\d2\\u', title => $anode, xrange => [0,255], yrange => [0,255], symbol => -1 };


}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
