package Network;

use POSIX;

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my %params = @_;
  my $self  = {};
  $self->{'file'}  = undef;
  $self->{'cachesize'}  = 1000;
  $self->{'blocksize'}  = 5000;
  $self->{'mode'}  = "r";

  $self->{'list'}  = {};
  $self->{'fh'}  = undef;
  $self->{'cache'}  = {};
  $self->{'num'}  = 0;
  $self->{'numBits'}  = 0;
  $self->{'numBytes'}  = 0;
  $self->{'startPtr'} = 0;

  foreach (keys(%params)) {
    $self->{$_} = $params{$_};
  }

  bless ($self, $class);
  return $self;
}

sub init {
  my $self = shift;
  my $file = $self->{'file'};
  if ($self->{'mode'} eq "r") {
    my $fh;
    open($fh, "<$file") || die "Can't open $file\n";
    binmode($fh);
    $self->{'fh'} = $fh;
  }
  if ($self->{'mode'} eq "rw") {
    my $fh;
    open($fh, "+<$file") || die "Can't open $file\n";
    binmode($fh);
    $self->{'fh'} = $fh;
  }
  if ($self->{'mode'} eq "w") {
    my $fh;
    open($fh, ">$file") || die "Can't open $file\n";
    binmode($fh);
    $self->{'fh'} = $fh;
  }
}

sub readHeader {
  my $self = shift;
  my $fh = $self->{'fh'};
  read($fh, $buffer, 1); 
  my $magic = unpack("C", $buffer);
  $self->{'magic'} = $magic; 
  read($fh, $buffer, 1); 
  my $major = unpack("C", $buffer);
  $self->{'major'} = $major; 
  read($fh, $buffer, 1); 
  my $minor = unpack("C", $buffer);
  $self->{'minor'} = $minor; 
}

sub readList {
  my ($self, $num) = @_;
  my $fh = $self->{'fh'};
  my $buffer;
  seek($fh, 3 + $num * 4, 0);
  read($fh, $buffer, 4);
  my $ptr = unpack("I", $buffer);
  seek($fh, $ptr, 0);
  read($fh, $buffer, 4);
  my $len = unpack("I", $buffer);
  read($fh, $buffer, $len);
  my $name = $buffer;
  read($fh, $buffer, 4);
  my $size = unpack("I", $buffer);
  my $res = [];
  for (my $i =0; $i < $size; $i ++) {
    read($fh, $buffer, 4);
    my $len = unpack("I", $buffer);
    read($fh, $buffer, $len);
    my $n = $buffer;
    push @$res, $n;
  }
  $self->{'list'}->{$name} = $res;
  return $res;
}

sub getLow {
  my $self = shift;
  return $self->readList(0);
}

sub getHigh {
  my $self = shift;
  return $self->readList(1);
}

sub getBalanced {
  my $self = shift;
  return $self->readList(2);
}

sub getNum {
  my $self = shift;
  return $self->{'num'};
}

sub getNumBits {
  my $self = shift;
  return $self->{'numBits'};
}

sub getNumBytes {
  my $self = shift;
  return $self->{'numBytes'};
}

sub getStartPtr {
  my $self = shift;
  return $self->{'startPtr'};
}

sub getMatrixEnd {
  my $self = shift;
  return $self->{'startPtr'} + $self->{'num'} * $self->{'numBytes'};
}

sub setMatrixEnd {
  my $self = shift;
  my $ptr = $self->getMatrixEnd();
  my $fh = $self->{'fh'};
  seek($fh, $ptr, 0);
}

sub readNetwork {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $buffer;

  seek($fh, 3 + 3 * 4, 0);
  read($fh, $buffer, 4);
  my $ptr = unpack("I", $buffer);
  read($fh, $buffer, 4);
  my $num = unpack("I", $buffer);
  read($fh, $buffer, 4);
  my $numBits = unpack("I", $buffer);
  my $numBytes = floor($num * $numBits/8) + 1;
  seek($fh, $ptr, 0);
  $self->{'num'} = $num;
  $self->{'numBits'} = $numBits;
  $self->{'numBytes'} = $numBytes;
  $self->{'startPtr'} = $ptr;
  #print "$num, $numBits, $numBytes, $ptr\n";
}

sub readNetworkFile {
  my $self = shift;
  $self->init();
  $self->readHeader();
  $self->readList(0);
  $self->readList(1);
  $self->readList(2);
  $self->readNetwork();
}

sub hashCode {
  my ($self, $a, $b) = @_;
  my $hash = $a * $self->{'numBytes'} + floor($b * $self->{'numBits'}/8);
  return floor($hash/$self->{'blocksize'});
}

sub loadBlock {
  my ($self, $hash)= @_;
  if (scalar(keys(%{$self->{'cache'}})) >= $self->{'cachesize'}) {
    $self->flushCache();
  }
  if (!defined $self->{'cache'}->{$hash}) {
    print STDERR "L $hash\n";
    my $pos = $self->{'startPtr'} + $hash * $self->{'blocksize'};
    my $buffer;
    my $fh = $self->{'fh'};
    seek($fh, $pos, 0);
    read($fh, $buffer, $self->{'blocksize'});
    $self->{'cache'}->{$hash} = $buffer;
  }
}

sub flushCache {
  my $self = shift;
  if ($self->{'mode'} eq "w" || $self->{'mode'} eq "rw") {
    foreach (sort keys(%{$self->{'cache'}})) {
      my $fh = $self->{'fh'};
      my $buffer = $self->{'cache'}->{$_};
      my $pos = $self->{'startPtr'} + $_ * $self->{'blocksize'};
      if ($pos != tell($fh)) {
        seek($fh, $pos, 0);
      }
      print $fh $buffer;
    }
    print STDERR "Cache clear\n";
  }
  $self->{'cache'} = {};
}

sub close {
  my $self = shift;
  $self->flushCache();
  my $fh = $self->{'fh'};
  close($fh);
}

sub writeHeader {
  my $self = shift;
  my $fh = $self->{'fh'};
  $self->{'magic'} = 0x55;
  $self->{'major'} = 1;
  $self->{'minor'} = 0;
  print $fh pack("C", $self->{'magic'});
  print $fh pack("C", $self->{'major'});
  print $fh pack("C", $self->{'minor'});
  for (my $i = 0; $i < 10; $i++) {
    $self->writeInt(0);
  }
}

sub writeHeader_1_1 {
  my $self = shift;
  my $fh = $self->{'fh'};
  $self->{'magic'} = 0x55;
  $self->{'major'} = 1;
  $self->{'minor'} = 1;
  print $fh pack("C", $self->{'magic'});
  print $fh pack("C", $self->{'major'});
  print $fh pack("C", $self->{'minor'});
  for (my $i = 0; $i < 10; $i++) {
    $self->writeInt(0);
  }
}

sub writeInt {
  my ($self, $val) = @_;
  my $fh = $self->{'fh'};
  my $str = pack("L", $val);
  print $fh  $str;
}

sub readInt {
  my $self = shift;
  my $fh = $self->{'fh'};
  my $buffer;
  read($fh, $buffer, 4);
  my $val = unpack("I", $buffer);
  return $val;
}

sub writeLengthPrefixString {
  my ($self, $str) = @_;
  my $fh = $self->{'fh'};
  $self->writeInt(length($str));
  print $fh  $str;
}

sub writeList {
  my ($self, $num, $tag, $list) = @_;
  my $fh = $self->{'fh'};
  my $ptr = tell($fh);
  seek($fh, 3 + $num * 4, 0);
  $self->writeInt($ptr);
  seek($fh, $ptr, 0);
  $self->writeLengthPrefixString($tag);
  my $size = scalar(@{$list});
  $self->writeInt($size);
  for (my $i =0; $i < $size; $i ++) {
    $self->writeLengthPrefixString($list->[$i]);
  }
  $self->{'list'}->{$tag} = $list;
}

sub startMatrix {
  my ($self, $num, $numBits) = @_;
  my $numBytes = floor($num * $numBits/8) + 1;
  my $fh = $self->{'fh'};
  $self->{'num'} = $num;
  $self->{'numBits'} = $numBits;
  $self->{'numBytes'} = $numBytes;
  $self->{'startPtr'} = tell($fh);
  seek($fh, 3 + 3 * 4, 0);
  $self->writeInt($self->{'startPtr'});
  $self->writeInt($self->{'num'});
  $self->writeInt($self->{'numBits'});
  seek($fh, $self->{'startPtr'}, 0);
}

sub writeCode {
  my ($self, $a, $b, $code) = @_;
  my $hash = $self->hashCode($a, $b);
  my $buffer;
  $self->loadBlock($hash);
  $buffer = $self->{'cache'}->{$hash};
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  my $byte_offset = ($a * $numBytes + floor($b * $numBits/8)) % $blocksize;
  my $bit_offset = ($b * $numBits) % 8;
  my $mask = (1 << $numBits) - 1;
  if (($bit_offset + $numBits) > 8) {
    if (($byte_offset+1) < $blocksize) {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      my $val1 = unpack "c", substr($buffer, $byte_offset+1, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $val = (($val & ~($mask << $bit_offset)) | (($code & $mask) << $bit_offset));
      substr($buffer, $byte_offset, 1) = pack "c", ($val & 0xff);
      substr($buffer, $byte_offset+1, 1) = pack "c", ($val >> 8);
    }
    else {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      $self->loadBlock($hash+1);
      $buffer = $self->{'cache'}->{$hash+1};
      my $val1 = unpack "c", substr($buffer, 0, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $val = (($val & ~($mask << $bit_offset)) | (($code & $mask) << $bit_offset));
      substr($buffer, 0, 1) = pack "c", ($val >> 8);
      $self->{'cache'}->{$hash+1} = $buffer;
      $self->loadBlock($hash);
      $buffer = $self->{'cache'}->{$hash};
      substr($buffer, $byte_offset, 1) = pack "c", ($val & 0xff);
    }
  }
  else {
    my $val = unpack "c", substr($buffer, $byte_offset, 1);
    $val = $val & 0xff;
    $val = (($val & ~($mask << $bit_offset)) | (($code & $mask) << $bit_offset));
    substr($buffer, $byte_offset, 1) = pack "c", ($val & 0xff);
    #my $pos = $startPtr + $hash * $blocksize;
    #print "$val $code $bit_offset $byte_offset $pos\n";
  }
  $self->{'cache'}->{$hash} = $buffer;
}

sub readCode {
  my ($self, $a, $b) = @_;
  my $hash = $self->hashCode($a, $b);
  my $buffer;
  $self->loadBlock($hash);
  $buffer = $self->{'cache'}->{$hash};
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  my $byte_offset = ($a * $numBytes + floor($b * $numBits/8)) % $blocksize;
  my $bit_offset = ($b * $numBits) % 8;
  my $mask = (1 << $numBits) - 1;
  my $code = 0;
  if (($bit_offset + $numBits) > 8) {
    if (($byte_offset+1) < $blocksize) {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      my $val1 = unpack "c", substr($buffer, $byte_offset+1, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $code = ($val >> $bit_offset) & $mask;
    }
    else {
      my $val0 = unpack "c", substr($buffer, $byte_offset, 1);
      $val0 = $val0 & 0xff;
      $self->loadBlock($hash+1);
      $buffer = $self->{'cache'}->{$hash+1};
      my $val1 = unpack "c", substr($buffer, 0, 1);
      $val1 = $val1 & 0xff;
      my $val = $val0 | ($val1 << 8);
      $code = ($val >> $bit_offset) & $mask;
    }
  }
  else {
    my $val = unpack "c", substr($buffer, $byte_offset, 1);
    $val = $val & 0xff;
    $code = ($val >> $bit_offset) & $mask;
    #my $pos = $startPtr + $hash * $blocksize;
    #print "$val $code $bit_offset $byte_offset $pos\n";
  }
  return $code;
}

sub print_1_0 {
  my $self = shift;
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  foreach my $name ("low", "high", "balanced") {
    my $size = scalar(@{$self->{'list'}->{$name}});
    print "$name ($size):";
    for (my $i =0; $i < $size; $i ++) {
      if (($i % 5) == 0) {
        print "\n";
      }
      my $id = $self->{'list'}->{$name}->[$i];
      print $id, ", ";
    }
    print "\n";
  }
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  my $balanced = $self->{'list'}->{"balanced"};
  for (my $i = 0; $i < $self->{'num'}; $i++) {
    print STDERR $i, "\n";
    for (my $j = 0; $j < $self->{'num'}; $j++) {
      my $code = $self->readCode($i, $j);
      if ($code > 0) {
        print "$code\t$i\t$j\t".$balanced->[$i]."\t".$balanced->[$j]."\n";
      }
    }
  }
}

# 0 - No relation
# 1 - $i low -> $j high
# 2 - $i low -> $j low
# 3 - $i high -> $j high
# 4 - $i high -> $j low
# 5 - Equivalent
# 6 - Opposite
sub readCode_1_1 {
  my ($self, $i, $j) = @_;
  my $b0 = $self->readCode(2 * $i,   2 * $j);
  my $b1 = $self->readCode(2 * $i,   2 * $j+1);
  my $b2 = $self->readCode(2 * $i+1, 2 * $j);
  my $b3 = $self->readCode(2 * $i+1, 2 * $j+1);
  my $total = $b0 + $b1 + $b2 + $b3;
  if ($total == 1) {
    if ($b0 == 1) { return 1; }
    if ($b1 == 1) { return 2; }
    if ($b2 == 1) { return 3; }
    if ($b3 == 1) { return 4; }
  }
  if ($total == 2) {
    if ($b1 == 1 && $b2 == 1) { return 5; }
    if ($b0 == 1 && $b3 == 1) { return 6; }
  }
  return 0;
}

sub hasStats {
  my $self = shift;
  my $buffer;
  my $fh = $self->{'fh'};
  seek($fh, 3 + 6 * 4, 0);
  read($fh, $buffer, 4);
  my $ptr = unpack("I", $buffer);
  return $ptr;
}

sub printStats {
  my $self = shift;
  my $balanced = $self->{'list'}->{"balanced"};
  print "Stats : ", $self->getMatrixEnd(), "\n";
  $self->setMatrixEnd();
  for (my $i = 0; $i < $self->{'num'}/2; $i++) {
    print $balanced->[$i];
    for (my $j = 0; $j < 7; $j++) {
      print "\t".$self->readInt();
    }
    for (my $j = 0; $j < 3; $j++) {
      $self->readInt();
    }
    print "\n";
  }
}

sub print_1_1 {
  my $self = shift;
  my $numBytes = $self->{'numBytes'};
  my $numBits = $self->{'numBits'};
  my $num = $self->{'num'};
  my $blocksize = $self->{'blocksize'};
  my $cachesize = $self->{'cachesize'};
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  foreach my $name ("low", "high", "balanced") {
    my $size = scalar(@{$self->{'list'}->{$name}});
    print "$name ($size):";
    for (my $i =0; $i < $size; $i ++) {
      if (($i % 5) == 0) {
        print "\n";
      }
      my $id = $self->{'list'}->{$name}->[$i];
      print $id, ", ";
    }
    print "\n";
  }
  print "Matrix [ n=$num, bit=$numBits, row=$numBytes bytes]\n";
  my $balanced = $self->{'list'}->{"balanced"};
  for (my $i = 0; $i < $self->{'num'}/2; $i++) {
    print STDERR $i, "\n";
    for (my $j = 0; $j < $self->{'num'}/2; $j++) {
      my $code = $self->readCode_1_1($i, $j);
      if ($code > 0) {
        print "$code\t$i\t$j\t".$balanced->[$i]."\t".$balanced->[$j]."\n";
      }
    }
  }
  if ($self->hasStats() != 0) {
    $self->printStats();
  }
}

sub print {
  my $self = shift;
  my $version = $self->{'major'}.".".$self->{'minor'};
  print "Version : $version\n";
  if ($version eq "1.0") {
    return $self->print_1_0();
  }
  if ($version eq "1.1") {
    return $self->print_1_1();
  }
}

sub readMatrix {
  my $self = shift;
  $self->readNetwork();
  my $fh = $self->{'fh'};
  my $numBytes = $self->{'numBytes'};
  my $num = $self->{'num'};
  my $network = [];
  my $buffer;
  for (my $i = 0; $i < $num; $i++) {
    print STDERR $i, "\n";
    read($fh, $buffer, $numBytes);
    if (length($buffer) < $numBytes) {
      my $str = "\0" x ($numBytes - length($buffer));
      $buffer = $buffer.$str;
    }
    $network->[$i] = $buffer;
  }
  return $network;
}

sub readIndex {
  my $indexFile = shift;
  my $idhash = {} ;
  my $namehash = {} ;
  my $revidhash = {} ;
  my $revnamehash = {} ;
  my $revptrhash = {} ;
  open(FL, "<$indexFile") || die "Can't open $indexFile\n";
  <FL>; # Header
  my $index = 0;
  while (<FL>) {
    my ($id, $name, $ptr, $desc) = split("\t", $_);
    my $idu = uc($id);
    my $nameu = uc($name);
    if (!defined $idhash->{$idu}) {
      $idhash->{$idu} = [];
    }
    push @{$idhash->{$idu}}, $index;
    if (!defined $namehash->{$nameu}) {
      $namehash->{$nameu} = [];
    }
    push @{$namehash->{$nameu}}, $index;
    $revptrhash->{$index} = $ptr;
    $revidhash->{$index} = $id;
    $revnamehash->{$index} = $name;
    $index++;
  }
  close(FL);
  $indexFile =~ s/.idx$/.thr/g;
  #my $thrHash = &readThrFile($indexFile);
  my $thrHash = {};
  $indexFile =~ s/.thr$/.pcl/g;
  return [$idhash, $namehash, $revidhash, $revnamehash, $revptrhash, $thrHash,
         $indexFile];
}

sub getLink {
  my ($gene, $org) = @_;
  if ($org eq "Mm") {
    return "http://smd.stanford.edu/cgi-bin/source/sourceResult?criteria=$gene&choice=Gene&option=Name&organism=Mm";
  }
  elsif ($org eq "Dm") {
    return "http://flybase.org/cgi-bin/uniq.html?context=$gene&species=Dmel&db=fbgn&caller=quicksearch";
  }
  elsif ($org eq "Sgd") {
    return "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=cerevisiae&name=$gene&desc=yes&submit=Search";
  }
  elsif ($org eq "Ath") {
    return "http://www.arabidopsis.org/servlets/TairObject?name=$gene&type=gene";
  }
  elsif ($org eq "Affy") {
    return "https://www.affymetrix.com/LinkServlet?probeset=$gene";
  }
  else {
    return "http://smd.stanford.edu/cgi-bin/source/sourceResult?criteria=$gene&choice=Gene&option=Name&organism=Hs";
  }
}

sub getSymbol {
  my ($index, $hashTable, $balanced) = @_;
  my $low = 0;
  if ($index =~ /b/) {
    $low = 1;
  }
  $index =~ s/b//g;
  my $probeid = $balanced->[$index];
  my $n = $hashTable->[0]->{uc($probeid)}->[0];
  my $name = $hashTable->[3]->{$n};
  if ($low == 0) {
    return "$name\-$probeid";
  }
  else {
    return "$name\_b-$probeid";
  }
}

sub readThrFile {
  my $thrFile = shift;
  my $hash = {};
  open(FL, "< $thrFile") || die "Can't open $thrFile\n";
  while (<FL>) {
    s/[\r\n]//g;
    my ($index, $thr1, $stat, $thr0, $thr2, @rest) = split("\t");
    $hash->{$index} = [$thr1, $stat, $thr0, $thr2];
  }
  return $hash;
}

sub getSplitName {
  my ($name, $len) = @_;
  my $cname = $name;
  if (length($name) > $len) {
    $cname = "";
    for ($i = 0; $i < length($name); $i+=$len) {
      $cname = $cname . substr($name, $i, $len) . "<br/>";
    }
  }
  return $cname;
}

sub printSubgraph {
  my ($self, $indexFile, $list, $org, $outdir) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version eq "1.0") {
    print "Only Relation format 1.1 is supported\n";
    exit;
  }
  print "Creating directory ...\n";
  if (! -e $outdir) {
    mkdir($outdir) || die "Can't create $outdir\n";
  }
  print "Reading Symbols ...\n";
  my $hashTable = &readIndex($indexFile);
  print "Done. Time=".times()."\n";
  print "Reading list ...\n";
  my $genes = {};
  my $fh;
  open($fh, "<$list") || die "Can't open $list\n";
  while (<$fh>) {
    s/[\r\n]//g;
    my @list = split(/[^\w-]+/);
    foreach my $w (@list) {
      $w = uc($w);
      next if ($w =~ /^-*$/);
      next if ($w eq "UNKNOWN");
      if (defined $hashTable->[0]->{$w}) {
        foreach my $n (@{$hashTable->[0]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
      if (defined $hashTable->[1]->{$w}) {
        foreach my $n (@{$hashTable->[1]->{$w}}) {
          $genes->{$n} = 1;
        }
      }
    }
  }
  close($fh);
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $lowgenes = {};
  my $highgenes = {};
  my $badgenes = {};
  my $listHash = {};
  my ($numbal, $numlow, $numhigh, $numbad) = (0, 0, 0, 0);
  foreach my $k (keys(%{$genes})) {
    my $probeid = $hashTable->[2]->{$k};
    my $name = $hashTable->[3]->{$k};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    if (defined $lowhash->{$probeid}) {
      push @{$lowgenes->{$name}}, $k;
      $numlow++;
    }
    elsif (defined $highhash->{$probeid}) {
      push @{$highgenes->{$name}}, $k;
      $numhigh++
    }
    elsif (defined $balancedhash->{$probeid}) {
      push @{$listHash->{$name}}, $k;
      $numbal++;
    }
    else {
      push @{$badgenes->{$name}}, $k;
      $numbad++;
    }
  }
  print "Done. Time=".times()."\n";
  print "Creating HTML files ... \n";
  my $indexfh;
  my $genesfh;
  my $gmlfh;
  open($indexfh, ">$outdir/index.html") || die "Can't open $outdir/index.html\n";
  open($genesfh, ">$outdir/genes.txt") || die "Can't open $outdir/genes.txt\n";
  open($gmlfh, ">$outdir/graph.gml") || die "Can't open $outdir/graph.gml\n";
  print $indexfh "<html> <body> <center>\n";
  print $indexfh "<h1>  Boolean interactions of genes </h1>\n";
  print $indexfh "<a href=\"genes.txt\"> List of genes </a> <br/>\n";
  print $indexfh "<a href=\"graph.gml\"> Graph in GML format</a> <br/>\n";

  print $indexfh "<table border=0>";
  print $indexfh "<tr> <td> Affy ID </td>\n";
  print $indexfh "<td> Total </td> <td> lohi </td>\n";
  print $indexfh "<td> lolo </td> <td> hihi </td> <td> hilo </td>\n";
  print $indexfh "<td> equ </td> <td> opo </td> <td> link </td> <td>Symbol</td></tr>\n";
  print $genesfh "balanced\t$numbal\n";
  print $genesfh "Affy ID\tTotal\tlohi\tlolo\thihi\thilo\tequ\topo\tSymbol\n";
  foreach my $name (sort keys(%{$listHash})) {
    foreach my $index (@{$listHash->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $relations = {};
      my $numrel = 0;
      my $i = $balancedhash->{$probeid};
  #print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
      foreach my $name1 (sort keys(%{$listHash})) {
        foreach my $index1 (@{$listHash->{$name1}}) {
          my $probeid1 = $hashTable->[2]->{$index1};
          my $j = $balancedhash->{$probeid1};
          my $code = $self->readCode_1_1($i, $j);
          if ($code > 0) {
            if (!defined $relations->{$code}) {
              $relations->{$code} = [];
            }
            push @{$relations->{$code}}, $j;
  #print "$code ".$balancedhash->{$id}." $i\n";
            $numrel++;
          }
        }
      }
      print "$probeid -> $numrel\n";
      next if ($numrel <= 0);
      my $cname = &getSplitName($name, 20);
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $namelink = &getLink($genename, $org);
      my $idlink = &getLink($probeid, "Affy");
      my $fileprobeid = $probeid;
      $fileprobeid =~ s/\//_/g;
      $fileprobeid =~ s/\:/_/g;
      print $genesfh "$probeid\t";
      print $genesfh "$numrel\t";
      print $genesfh scalar(@{$relations->{1}})+0,"\t";
      print $genesfh scalar(@{$relations->{2}})+0,"\t";
      print $genesfh scalar(@{$relations->{3}})+0,"\t";
      print $genesfh scalar(@{$relations->{4}})+0,"\t";
      print $genesfh scalar(@{$relations->{5}})+0,"\t";
      print $genesfh scalar(@{$relations->{6}})+0,"\t";
      print $genesfh "$name\n";
      print $indexfh "<tr>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$idlink\"> $probeid </a> </td>";
      print $indexfh "<td> $numrel </td>";
      print $indexfh "<td>", scalar(@{$relations->{1}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{2}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{3}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{4}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{5}})+0,"</td>";
      print $indexfh "<td>", scalar(@{$relations->{6}})+0,"</td>";
      print $indexfh "<td> <a href=\"$fileprobeid\.html\"> link </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$namelink\"> $cname </a> </td>";
      print $indexfh "</tr>";
      my $idfh;
      open($idfh, ">$outdir/$fileprobeid\.html") || die "Can't open $outdir/$fileprobeid\.html\n";
print $idfh <<END;
<html>
<head>
<SCRIPT SRC="http://gourd.stanford.edu/~sahoo/public/mktree.js"
LANGUAGE="JavaScript"></SCRIPT>
<LINK REL="stylesheet"
HREF="http://gourd.stanford.edu/~sahoo/public/mktree.css">
</head>
<body>
<center>
<h1> Boolean interactions of genes </h1>
<br>
<A href="#" onClick="expandTree('tree1'); return
false;">ExpandAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="collapseTree('tree1'); return
false;">CollapseAll</A>&nbsp;&nbsp;&nbsp;
<A href="#" onClick="expandTreeDepth('tree1',2); return false;">ExpandDepth 2
</A>

</center>
<ul class="mktree" id="tree1">
END
      print $idfh "<li> <font size=+2> <a target=\"_blank\" href=\"$namelink\"> $genename </a> </font> - <a target=\"_blank\" href=\"$idlink\"> $probeid </a> (", $numrel ,")<ul>\n";
      &printListHtml($relations->{5}, $index, "Equivalent", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{6}, $index, "Opposite", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{2}, $index, "($genename low -> B low)", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{1}, $index, "($genename low -> B high)", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{4}, $index, "($genename high -> B low)", $org, $hashTable, $balanced, $idfh);
      &printListHtml($relations->{3}, $index, "($genename high -> B high)", $org, $hashTable, $balanced, $idfh);
      print $idfh "</ul> </body> </html>\n";
    }
  }
  print $indexfh "</table><br/>";

  print $indexfh "<b> High Genes </b>($numhigh)<br/>\n";
  print $indexfh "<table border=0>\n";
  print $genesfh "high\t$numhigh\n";
  foreach my $name (sort keys(%{$highgenes})) {
    foreach my $index (@{$highgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print $genesfh "$probeid\t$name\n";
      print $indexfh "<tr> <td> <a target=\"_blank\" href=\"$ilink\"> $probeid </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$link\"> $cname </a> </td> </tr>\n";
    }
  }
  print $indexfh "</table><br/>\n";
  print $indexfh "<b> Low Genes </b>($numlow)<br/>\n";
  print $indexfh "<table border=0>\n";
  print $genesfh "low\t$numlow\n";
  foreach my $name (sort keys(%{$lowgenes})) {
    foreach my $index (@{$lowgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print $genesfh "$probeid\t$name\n";
      print $indexfh "<tr> <td> <a target=\"_blank\" href=\"$ilink\"> $probeid </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$link\"> $cname </a> </td> </tr>\n";
    }
  }
  print $indexfh "</table><br/>\n";
  print $indexfh "<b> Bad Dynamic Range Genes </b>($numbad)\n";
  print $indexfh "<table border=0>\n";
  print $genesfh "bad\t$numbad\n";
  foreach my $name (sort keys(%{$badgenes})) {
    foreach my $index (@{$badgenes->{$name}}) {
      my $probeid = $hashTable->[2]->{$index};
      my $ilink = &getLink($probeid, "Affy");
      my $genename = $name;
      $genename =~ s/\/.*//g;
      my $link = &getLink($genename, $org);
      my $cname = &getSplitName($name, 20);
      print $genesfh "$probeid\t$name\n";
      print $indexfh "<tr> <td> <a target=\"_blank\" href=\"$ilink\"> $probeid </a> </td>";
      print $indexfh "<td> <a target=\"_blank\" href=\"$link\"> $cname </a> </td> </tr>\n";
    }
  }
  print $indexfh "</table><br/>\n";
  print $indexfh "</center></body></html>\n";
  close($indexfh);
  close($genesfh);
  close($gmlfh);
  print "Done. Time=".times()."\n";
}

sub printText {
  my($self, $indexFile, $probeid, $type) = @_;
  my $version = $self->{'major'}.".".$self->{'minor'};
  if ($version ne "1.1") {
    print "Only version 1.1 supported\n";
    exit(1);
  }
   
  $probeid = uc($probeid);
  my $file = $self->{'file'};
  my $hashTable = &readIndex($indexFile);
  if (!defined $hashTable->[0]->{$probeid}) {
    print "Cannot fine probe $probeid\n";
    return;
  }
  my $low = $self->getLow();
  my $high = $self->getHigh();
  my $balanced = $self->getBalanced();

  my $lowhash = {};
  my $highhash = {};
  my $balancedhash = {};
  for(my $i = 0; $i < scalar(@{$low}); $i++) {
    $lowhash->{$low->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$high}); $i++) {
    $highhash->{$high->[$i]} = $i;
  }
  for(my $i = 0; $i < scalar(@{$balanced}); $i++) {
    $balancedhash->{$balanced->[$i]} = $i;
  }
  my $index = $hashTable->[0]->{$probeid}->[0];
  my $name = $hashTable->[3]->{$index};
  my $id = $hashTable->[2]->{$index};
  my $cname = &getSplitName($name, 20);
  my $genename = $name;
  $genename =~ s/\/.*//g;
  if (defined $lowhash->{$id}) {
    print "$id\t$genename\t(Always Low)\n";
    return;
  }
  if (defined $highhash->{$id}) {
    print "$id\t$genename\t(Always High)\n";
    return;
  }
  if (!defined $balancedhash->{$id}) {
    print "$id\t$genename\t(Bad dynamic range)\n";
    return;
  }
  my $ptr = $self->getStartPtr();
  my $num = $self->getNum();
  my $numBits = $self->getNumBits();
  my $numBytes = $self->getNumBytes();
  my $relations = {};
  my $numrel = 0;
  my $i = $balancedhash->{$id};
  #print "$num $numBits $numBytes $ptr ".$balancedhash->{$id}."\n";
  for (my $j = 0; $j < $num/2; $j++) {
    my $code = $self->readCode_1_1($i, $j);
    if ($code > 0) {
      if (!defined $relations->{$code}) {
        $relations->{$code} = [];
      }
      push @{$relations->{$code}}, $j;
      #print "$code ".$balancedhash->{$id}." $i\n";
      $numrel++;
    }
  }
  my @status = ("No relation", "($genename low -> B high)",
  "($genename low -> B low)", "($genename high -> B high)",
  "($genename high -> B low)", "Equivalent", "Opposite");
  my $list = $relations->{$type};  
  print "$id\t$genename\t$numrel\t", scalar(@{$list}), "\t", $status[$type], "\n";
  for (my $i = 0; $i < scalar(@{$list}); $i++) {
    my $probeid = $balanced->[$list->[$i]];
    my $nid = $hashTable->[0]->{uc($probeid)}->[0];
    my $name = $hashTable->[3]->{$nid};
    if (!defined $name || $name eq "---" || $name eq "-" || $name eq "") {
      $name = "Unknown";
    }
    print "$probeid\t$name\n";
  }
}

package main;

if (scalar(@ARGV) <= 1) {
    print "Usage : extract.pl <cmd> <args>\n";
    print "<type> :\n";
    print "       0 - No relation\n";
    print "       1 - A low  -> B high\n";
    print "       2 - A low  -> B low\n";
    print "       3 - A high -> B high\n";
    print "       4 - A high -> B low\n";
    print "       5 - Equivalent\n";
    print "       6 - Opposite\n";
    print "<cmd> :\n";
    print "   text <network.rl> <network.idx> <probeid> <type>\n";
    print "   index <pclfile>\n";
    exit;
}

use POSIX;

my $cmd = shift @ARGV;

if ($cmd eq "text") {
  my ($file, $indexFile, $probeid, $type) = @ARGV;
  my $n = Network->new(file => $file, mode => "r");
  $n->readNetworkFile(); # read headers
  $n->printText($indexFile, $probeid, $type);
}

if ($cmd eq "index") {
  &buildIndex(@ARGV);
}

sub buildIndex {
  my $pclFile = shift;

  print "ID\tName\tPtr\tDesc\n";
  open(FL, "< $pclFile") || die "Can't open $pclFile\n";
  my $index = 0;
  my $ptr = tell(FL);
  while (<FL>) {
    my @list = split("\t", $_);
    my ($name, $desc) = split(":", $list[1], 2);
    if ($name =~ /^\s*$/) {
       $name = "---";
    }
    if ($name eq "MGI") {
      my ($id, $rest) = split(":", $desc);
      $name = $name.":".$id;
      $desc = $rest;
    }
    print $list[0]."\t".$name."\t".$ptr."\t".$desc."\n";
    $index++;
    $ptr = tell(FL);
  }
  close(FL);
}

