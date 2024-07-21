package Lab;
use strict;

use PDL;
use PDL::Image2D;

require Exporter;
use vars qw( @ISA @EXPORT @EXPORT_OK $BGDIR $EVTDIR $TESTFILE $ANALDIR $TYPE %NBINS %BGFILES );

use MyRDB qw( rdb_cols );

@ISA = qw( Exporter );
@EXPORT_OK = qw(
		zero_hists
		samp_calc
		test_data
		rawy_limits_mcp
		rawy_limits_chipid
		CHIPX_MIN
		CHIPX_MAX
		CHIPY_MIN
		CHIPY_MAX
		RAWX_MIN
		RAWX_MAX
		RAWY_MIN
		RAWY_MAX
		TAPSIZE
		SUBTAPS
	       );

$BGDIR = '/data/legs/rpete/data/hrcs_lab/bg';
$EVTDIR = '/data/legs/rpete/data/hrcs_lab/evt1';
$ANALDIR = '/data/legs/rpete/data/hrcs_lab/analysis';
$TESTFILE = '/data/legs/rpete/cal/hrcs_gain/hrcs_lab.rdb';

# which binfile type to read/write
# samp, spimean, spimed or pha
$TYPE = 'samp';

# number of bins in the various histogram types
%NBINS = (
	  samp => 512,
	  spimean => 512,
	  spimed => 512,
	  pha => 256,
	 );

%BGFILES = (
	    samp => $ANALDIR.'/merged_bg_samp.bin',
	    pha => $ANALDIR.'/merged_bg_pha.bin',
	    spimean => $ANALDIR.'/merged_bg_spimean.bin',
	    spimed => $ANALDIR.'/merged_bg_spimed.bin',
	   );

use constant CHIPX_MIN => 1;
use constant CHIPX_MAX => 4096;

use constant CHIPY_MIN => 1;
use constant CHIPY_MAX => 16456;

use constant RAWX_MIN => 0;
use constant RAWX_MAX => 4095;

use constant RAWY_MIN => 0;
use constant RAWY_MAX => 16384 * 3 - 1;

use constant TAPSIZE => 256;
use constant SUBTAPS => 3;

sub merged_bg_omit {
  return qw(
	    p197061504
	    p197061508
	    p197061512
	    p197061516
	    );
}

sub samp_calc {
  my ($sumamps, $amp_sf) = @_;
  return $sumamps*(2**($amp_sf-1))/128;
#  return( $sumamps >> (8 - $amp_sf) );
}

# return all lines/energies in the lab configuration file, ordered by energy
sub lines {
  my ($lines, $energies) = test_data();

  my %lines; @lines{@$lines} = @$energies;

  my @lines = keys %lines;
  my @energies = values %lines;

  my @i = sort { $energies[$a] <=> $energies[$b] } 0..$#energies;

  return [@lines[@i]], [@energies[@i]];
}

# return energies for the given lines in the lab setup
sub energies {
  my @lines = @_;

  my ($l, $e) = lines();

  my %lines; @lines{@$l} = @$e;

  return map { $lines{$_} or die "line $_ not found"; $lines{$_} } @lines;
}

sub mcp_to_chipid {
  for ($_[0]) {
    $_ == 1  and return 1;
    $_ == 0  and return 2;
    $_ == -1 and return 3;
  }
  die "unrecognized MCP == $_[0]";
}

sub chipid_to_mcp {
  for ($_[0]) {
    $_ == 1 and return 1;
    $_ == 2 and return 0;
    $_ == 3 and return -1;
  }
  die "unrecognized chip_id == $_[0]";
}

sub rawy_limits_mcp {
  return rawy_limits_chipid(mcp_to_chipid(@_));
}

sub rawy_limits_chipid {
  my $id = shift;
  $id==1 or $id==2 or $id==3 or die "invalid chip_id == $id";
  return ((RAWY_MAX+1)/3 * ($id - 1), (RAWY_MAX+1)/3 * $id - 1);
}

sub chipy_limits_chipid {
}

sub test_data_hrcfile {
  my @hrcfiles = @_;

  # retrieve all test data, filter requested hrc file basenames
  my @keys = qw( line energy MCP time HRC_file b_time b_HRC_file );

  my %data;
  @data{@keys} = test_data();

  for my $f (@hrcfiles) {
    my @i = grep $data{HRC_file}[$_] eq $f, 0..$#{$data{HRC_file}};
    @i == 0 and die "no data found for basename $f";
    @i == 1 or die "multiple data found for basename $f";
    push @{$data{$_.'_retval'}}, @{$data{$_}}[@i] for @keys;
  }

  return map { $data{$_.'_retval'} } @keys;
}

sub test_data {

  my @anodes = @_;

  my ($line, $energy, $MCP, $time, $HRC_file, $b_time, $b_HRC_file) =
    rdb_cols($TESTFILE, qw( line energy MCP time HRC_file b_time b_HRC_file ) )
      or die;

  # extract entries matching specific anodes, if requested
  if (@anodes) {
    my @i;
    my %done;
    for my $anode (@anodes) {
      $done{$anode} and die "duplicate line = $anode";
      $done{$anode} = 1;
      my @ii = grep { $line->[$_] eq $anode } 0..$#{$line}
	or die "no lines matching $anode";
      push @i, @ii;
    }
    @$_ = @{$_}[@i]
      for $line, $energy, $MCP, $time, $HRC_file, $b_time, $b_HRC_file;
  }

  return $line, $energy, $MCP, $time, $HRC_file, $b_time, $b_HRC_file;
}

sub merged_bg_exptimes {
  my ($MCP, $time, $HRC_file, $b_time) = (test_data())[2..5];
  my %times = (1 => 0, 2 => 0, 3 => 0);
  my @omit = merged_bg_omit();
  my %omit; @omit{@omit} = ();
  for my $i (0..$#{$MCP}) {
    my $data_chipid = mcp_to_chipid($MCP->[$i]);
    for my $chipid (keys %times) {
      $times{$chipid} += $b_time->[$i];
      if (!exists $omit{$HRC_file->[$i]}) {
	$times{$chipid} += $time->[$i] unless $chipid == $data_chipid;
      }
    }
  }
  return %times;
}

sub zero_hists {
  my ($nx, $ny) = dims2d();

  my $hists = zeroes(float, $NBINS{$TYPE}, $nx, $ny);
  return $hists;
}

=begin comment

chip	crsu	crsv		no iffy crsv
----	----	----		------------
1	4b:13b	6b:63c		6b:63c
2	4c:13b	66b:125c	66b:125b
3	5b:13b	128b:186a	129a:186a

BG	2b:13b	see above	see above

=cut

BEGIN {

  my @usrc = qw( 4b:13b 4c:13b 5b:13b );
  my @ubg = ('2b:13b')x@usrc;
  my @viffy =   qw( 6b:63c 66b:125c 128b:186a );
  my @vnoiffy = qw( 6b:63c 66b:125b 129a:186a );

  sub active_mask_src {
    return active_mask_uv( \(@usrc, @viffy) );
  }
  sub active_mask_src_noiffyv {
    return active_mask_uv( \(@usrc, @vnoiffy) );
  }
  sub active_mask_bg {
    return active_mask_uv( \(@ubg, @viffy) );
  }
  sub active_mask_bg_noiffyv {
    return active_mask_uv( \(@ubg, @vnoiffy) );
  }
};

sub dims2d {
  my $nx = (RAWX_MAX - RAWX_MIN + 1) * SUBTAPS / TAPSIZE;
  my $ny = (RAWY_MAX - RAWY_MIN + 1) * SUBTAPS / TAPSIZE;
  return $nx, $ny;
}

sub active_mask_uv {
  my ($u, $v) = @_;
  my ($nx, $ny) = dims2d();

  my $mask = zeroes(byte, $nx, $ny);

  for my $i (0..$#{$u}) {
    my $u = $u->[$i];
    my $v = $v->[$i];

    (my $tmp = $mask->mslice([bradto2d($u)], [bradto2d($v)])) .= 1;
  }
  return $mask;
}

sub bradto2d {
  return map {
    my ($tap, $sub) = /(\d+)([abc])/ or die;
    $tap * SUBTAPS + ord($sub) - ord('a');
  } split ':', $_[0];
}

sub three_by_three {
  my $hists = shift;
  return $hists->xchg(0,2)->conv2d(ones($hists->type,3,3),{Boundary=>'Truncate'})->xchg(0,2)->sever;
}


sub bg_hists {
  my $file = @_ ? shift : $BGFILES{$TYPE};

  my $hists = zero_hists();
  add_to_hists($hists, $file);

  return $hists;
}

sub src_hists {
  my $anode = shift;
  my $yfilter = @_ ? shift : 1; # filter on rawy by default

  my $hists = zero_hists();

  my ($MCP, $HRC_file) = (Lab::test_data($anode))[2,4];

  for (0..$#{$HRC_file}) {
    my $hrcfile = $HRC_file->[$_];
    my $file = "$Lab::ANALDIR/${hrcfile}_$TYPE.bin";
    add_to_hists($hists, $file,
		 ($yfilter ? (rawy_limits_mcp($MCP->[$_])) : ()),
		 );
  }

  return $hists;
}

sub add_to_hists {
  my ($hists, $file) = (shift, shift);

  my ($yfilter, $rawy_lo, $rawy_hi);
  if (@_) {
    $yfilter = 1;
    ($rawy_lo, $rawy_hi) = (shift, shift);
  }

  my $it = Lab::InBinFile->new($file) or die;
  while (my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y) = $it->next_subtap) {
    next if ( $yfilter and ($y1 < $rawy_lo or $y2 > $rawy_hi) );
    my $xcoord = $xtap * SUBTAPS + $xsubtap;
    my $ycoord = $ytap * SUBTAPS + $ysubtap;
    (my $tmp = $hists->slice(",($xcoord),($ycoord)")) += $y;
  }
}

#
# images of PHA stats
#
# ideally we could just do this, but since rld() pads the first output
# dimension it's incorrect
#my $pha = $hists->rld(sequence(long,256))->float;
#my $stats = $pha->statsover;
#
# so brute force for now
#
sub hists_stats {
  my $hists = shift;

  my ($nz, $nx, $ny) = $hists->dims;

  my $xvals = sequence($hists->type, $nz);

  my $mean = zeroes(float, $nx, $ny);
  my $prms = $mean->copy;
  my $median = $mean->copy;
  my $min = $mean->copy;
  my $max = $mean->copy;
  my $adev = $mean->copy;
  my $rms = $mean->copy;

  for my $x (0..$nx-1) {
    for my $y (0..$ny-1) {
      my $hist = $hists->slice(",($x),($y)");
      next unless $hist->sum;
      my $vals = $hist->rld($xvals)->float;
      my @s = $vals->stats;
      $mean->set($x,$y,$s[0]);
      $prms->set($x,$y,$s[1]);
      $median->set($x,$y,$s[2]);
      $min->set($x,$y,$s[3]);
      $max->set($x,$y,$s[4]);
      $adev->set($x,$y,$s[5]);
      $rms->set($x,$y,$s[6]);
    }
  }

  return ($mean, $prms, $median, $min, $max, $adev, $rms);
}

sub hists_prob {
  my $hists = shift;
  my $n = $hists->getdim(0);
  return $hists / $hists->sumover->dummy(0);
}

sub hists_mean {
  my $hists = shift;
  my $xvals = @_ ? shift : sequence($hists->type, $hists->getdim(0));

  return ($hists * $xvals)->sumover / $hists->sumover;

  my $prob = hists_prob($hists);
  my $exp_x = ( $xvals * $prob )->sumover;
  return $exp_x;
}

sub hists_prms {
  my $hists = shift;
  my $xvals = @_ ? shift : sequence($hists->type, $hists->getdim(0));
  my $exp_x = hists_mean($hists, $xvals);
  my $prms = sqrt(
		  ($hists * ( $xvals - $exp_x->dummy(0) )**2 )->sumover
		    / ($hists->sumover-1)
		   );
  return $prms;
}

sub hists_rms {
  my ($hists, $xvals) = @_;

  my $exp_x2 = ( $xvals * $xvals * hists_prob($hists) )->sumover;
  my $exp_x = hists_mean($hists, $xvals);
  return sqrt($exp_x2 - $exp_x * $exp_x);
}

sub hists_median {
  my $hists = shift;
  my ($nz, $nx, $ny) = $hists->dims;

  my $xvals = sequence(long, $nz);

  my $median = zeroes(float, $nx, $ny);

  for my $x (0..$nx-1) {
    for my $y (0..$ny-1) {

      my $hist = $hists->slice(",($x),($y)");
      next unless $hist->sum;

      if (($hist->type == double) or ($hist->type == float)) {
	$median->set($x,$y,hist_median_float($xvals, $hist));
      }
      else {
	my $vals = $hist->rld($xvals)->float;
	$median->set($x,$y,$vals->median);
      }
    }
  }

  return $median;
}

# For background-subtracted data containing negative counts. Adapted
# from the hist_stats function in genstats.pl.
sub hist_median_float {
  my ($x, $y) = @_;

  my $y_sum = $y->sum;

  my $median = 0;
  eval {
    my $cumu_bounds = zeroes($x->nelem);
     (my $tmp = $cumu_bounds->slice('0:-2')) .= ($x->slice('0:-2') + $x->slice('1:-1'))/2;
    $cumu_bounds->set(-1, $x->at(-1)+($x->at(-1)-$x->at(-2))/2);
    my $cumu = $y->cumusumover / $y_sum;
    $median = interpol(0.5, $cumu, $cumu_bounds);
  };

  return $median;
}

# FIXME does not agree with hists_stats above since rld pads with zeroes
sub hists_min {
  my ($hists, $xvals) = @_;
  return $hists->rld($xvals)->minimum;
}

sub hists_max {
  my ($hists, $xvals) = @_;
  return $hists->rld($xvals)->maximum;
}

sub hists_adev {
  my ($hists, $xvals) = @_;

  my $n = $hists->getdim(0);

  my $prob = hists_prob($hists);
  my $exp_x = hists_mean($hists, $xvals);
  my $exp_adev = ( abs($xvals - $exp_x->dummy(0)) * $prob)->sumover;
  return $exp_adev;
}

sub quantile {
  @_ == 3 or @_ == 2 or die;
  my ($data, $f) = (shift, shift);
  my $sorted = @_ ? shift : 0;
  $data = $data->qsort unless $sorted;

  my $n = $data->getdim(0);
  my $i = rint(($n-1) * $f)->at;
  my $ii = $i+1;
  my $delta = ($n-1) * $f - $i;

  my $s1 = "($i)";
  my $s2 = "($ii)";
  return (1-$delta) * $data->slice($s1) + $delta * $data->slice($s2);
}

sub iqr {
  @_ == 1 or @_ == 2 or die;
  my $data = shift;
  my $sorted = @_ ? shift : 0;
  $data = $data->qsort unless $sorted;
  return quantile($data, 0.75, 1) - quantile($data, 0.25, 1);
}

sub trim {
  @_ == 2 or @_ == 3 or die;
  my ($data, $f) = (shift, shift);
  my $sorted = @_ ? shift : 0;

  $data = $data->qsort unless $sorted;

  my $i = int($f * $data->getdim(0));
  return $data->slice($i.':'.(-$i-1));
}

sub hists_trim {
  my $hists = (shift)->copy;
  my $f = shift;

  my $frac = ($hists->cumusumover / $hists->sumover->dummy(0))->flat;

  my $i = which(($frac<$f) | ($frac>1-$f));
  (my $tmp = $hists->flat->index($i)) .= 0;

  return $hists;
}

sub hists_indexND {
  my ($hists, $i) = @_;

  my $hdim = $hists->getdim(0);

  my $tmp;

  my $n = $i->getdim(1);

  my $ii = zeroes(long,3,$n);
  ($tmp = $ii->slice('1:-1,')) .= $i;
  $ii= $ii->dummy(1,$hdim)->clump(1,2)->copy;
  ($tmp=$ii->slice('0:0,,')->flat).= sequence(long,$hdim)->dummy(1,$n)->flat;

#  return $hists->indexND($ii)->reshape($hdim, $n); # makes copy
  return $hists->indexND($ii)->splitdim(0, $hdim);
}

package Lab::BinFile;
use strict;

use vars qw( @ISA @EXPORT_OK );

use Exporter;
@ISA = qw( Exporter );

use constant hist_vals => 1;
use constant data_vals => 2;

@EXPORT_OK = qw( hist_vals data_vals );

package Lab::InBinFile;
use strict;
use IO::File;
use PDL;

sub new {

  my $class = ref $_[0] ? ref shift : shift;

  my $filename = shift;

  my $type = @_ ? shift : $Lab::TYPE;

  my $fh = IO::File->new('< '.$filename) or die "error opening $filename: $!";

  my $ref = {
	     FH => $fh,
	     GOOD => 1,
	     TYPE => $type,
	     };

  bless $ref, $class;
}

sub next_subtap {
  my $ref = shift;
  my $fh = $ref->{FH};

  my $return_data = @_ ? shift : 1;

  my $hdr;
  $fh->read($hdr, 44) == 44 or $ref->{GOOD} = 0, return;

  my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2,
      $sample_type, $data_type, $data_n)
    = unpack('N*', $hdr);

  my ($type, $typesize);

  for ($data_type) {
    $_ == float(0)->get_datatype and $type=float(), $typesize=4, last;
    $_ == long(0)->get_datatype and $type=long(), $typesize=4, last;
    $_ == short(0)->get_datatype and $type=short(), $typesize=2, last;
    die $data_type;
  }

  my $hdim = $Lab::NBINS{$ref->{TYPE}};

  my ($x, $y);
  if ($return_data) {

    my $tmp;
    $fh->read($tmp, $data_n * $typesize) == $data_n * $typesize or die;

    my $data = zeroes($type, $data_n);
    ${ $data->get_dataref } = $tmp;

    $data->upd_data;
    $data->bswap4 unless isbigendian();

    for ($sample_type) {
      $_ == Lab::BinFile::hist_vals
	and $y = $data->long, last;
      $_ == Lab::BinFile::data_vals
	and $y = hist($data,-.5, $hdim-0.5,1)->long, last;
      die $sample_type
    };

    # e.g., SAMP histogram is in the file, but new() wasn't told
    die unless $hdim == $y->nelem;

    $x = sequence(long, $hdim);
  }

  # just advance the file pointer if we don't return the data
  else {
    $fh->seek($data_n * $typesize, SEEK_CUR) or die;
  }

  return $ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $x, $y;
}

package Lab::OutBinFile;
use strict;
use IO::File;
use PDL;

sub new {

  my $class = ref $_[0] ? ref shift : shift;

  my $filename = shift;
  my $type = @_ ? shift : $Lab::TYPE;

  my $fh = IO::File->new('> '.$filename) or die "error opening $filename for writing: $!";

  my $ref = {
	     FH => $fh,
	     TYPE => $type,
	     };

  bless $ref, $class;
}

sub add_subtap {
  my $ref = shift;
  my $fh = $ref->{FH};

  my ($ytap, $ysubtap, $xtap, $xsubtap, $y1, $y2, $x1, $x2, $vals) = @_;
  my $nbins = $Lab::NBINS{$ref->{TYPE}};

  my $data = $ref->{TYPE} ne 'pha' ? $vals->float->copy : $vals->short->copy;
  my $sample_type = Lab::BinFile::data_vals;

  if ($data->nelem > $nbins) {
    $data = hist($vals, -0.5, $nbins - 0.5, 1)->long;
    $sample_type = Lab::BinFile::hist_vals;
  }

  $data->bswap4 unless isbigendian();

  my $hdr = pack('N*',
		 $ytap, $ysubtap, $xtap, $xsubtap,
		 $y1, $y2, $x1, $x2,
		 $sample_type, $data->get_datatype, $data->nelem,
		 );

  $fh->print($hdr, ${ $data->get_dataref });
}

1;
