#!/usr/bin/perl -w

# Author: lh3

use strict;
use warnings;
use Getopt::Std;

my $version = '0.2.3';
&usage if (@ARGV < 1);

my $command = shift(@ARGV);
my %func = (snpFilter=>\&snpFilter, indelFilter=>\&indelFilter, showALEN=>\&showALEN, pileup2fq=>\&pileup2fq);

die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}};
exit(0);

#
# showALEN
#

sub showALEN {
  die(qq/Usage: samtools.pl showALEN <in.sam>\n/) if (@ARGV == 0 && -t STDIN);
  while (<>) {
	my @t = split;
	my $l = 0;
	$_ = $t[5];
	s/(\d+)[SMI]/$l+=$1/eg;
	print join("\t", @t[0..5]), "\t$l\t", join("\t", @t[6..$#t]), "\n";
  }
}

#
# indelFilter
#

sub indelFilter {
  my %opts = (D=>100, m=>10, r=>undef, s=>100);
  getopts('D:m:rs:', \%opts); # -s for scaling factor in score calculation

  die(qq/
Usage:   samtools.pl indelFilter [options] <in.indel>\n
Options: -D INT    maximum read depth [$opts{D}]
         -m INT    minimum distance between two adjacent indels [$opts{m}]
\n/) if (@ARGV == 0 && -t STDIN);

  my (@arr1, @arr2);
  my ($curr, $last) = (\@arr1, \@arr2);
  my $is_ref = defined($opts{r})? 1 : 0;
  while (<>) {
	my @t = split;
	next if ($t[2] ne '*');
	if (!$is_ref) {
	  next if ($t[3] eq '*/*');
	  next if ($t[5] == 0);
	}
	next if ($t[7] > $opts{D});
	# calculate indel score
	my $score = $t[5];
	$score += $opts{s} * $t[10] if ($t[8] ne '*');
	$score += $opts{s} * $t[11] if ($t[9] ne '*');
	@$curr = ($t[0], $t[1], $score, $_);
	my $do_swap = 1;
	if (defined $last->[0]) {
	  if ($curr->[0] eq $last->[0] && $last->[1] + $opts{m} > $curr->[1]) {
		$do_swap = 0 if ($last->[2] > $curr->[2]);
	  } else { # then print
		print $last->[3];
	  }
	}
	if ($do_swap) {
	  my $tmp = $curr; $curr = $last; $last = $tmp;
	}
  }
  print $last->[3] if (defined $last->[0]);
}

#
# snpFilter
#

sub snpFilter {
  my %opts = (f=>'', Q=>40, d=>3, w=>10, D=>0, N=>2, W=>10, q=>20, s=>50);
  getopts('f:s:w:q:Q:d:D:W:N:', \%opts);
  die(qq{
Usage:   samtools.pl snpFilter [options] <cns2snp.snp>

Options: -d INT        minimum depth to call a SNP [$opts{d}]
         -D INT        maximum depth, 0 to ignore [$opts{D}]
         -Q INT        required max mapping quality of the reads covering the SNP [$opts{Q}]
         -q INT        minimum SNP quality [$opts{q}]

         -f FILE       filtered samtools indels [null]
         -s INT        minimum samtols indel score [$opts{s}]
         -w INT        SNP within INT bp around an indel to be filtered [$opts{w}]

         -W INT        window size for filtering dense SNPs [$opts{W}]
         -N INT        maximum number of SNPs in a window [$opts{N}]
\n}) unless (@ARGV);
  my (%hash, $fh);
  my $skip = $opts{w};
  $opts{D} = 100000000 if ($opts{D} == 0);
  if ($opts{f}) { # filtered samtools indel
	my $n = 0;
	open($fh, $opts{f}) || die;
	while (<$fh>) {
	  my @t = split;
	  next if ($t[2] ne '*' || $t[3] eq '*/*' || $t[5] < $opts{s});
	  for (my $x = $t[1] - $skip + 1; $x < $t[1] + $skip; ++$x) {
		$hash{$t[0],$x} = 1;
	  }
	}
	close($fh);
  }
  my (@last, $last_chr);
  $last_chr = '';
  while (<>) {
	my @t = split;
	next if ($t[2] eq '*' || $hash{$t[0],$t[1]});
	next if ($t[2] eq $t[3]);
	my $is_good = ($t[7] >= $opts{d} && $t[7] <= $opts{D} && $t[6] >= $opts{Q} && $t[5] >= $opts{q})? 1 : 0;
	next unless ($is_good); # drop
	if ($t[0] ne $last_chr) { # a different chr, print
	  map { print $_->{L} if ($_->{F}) } @last;
	  @last = ();
	  $last_chr = $t[0];
	}
	# The following block implemented by Nathans Weeks.
	push(@last, {L => $_, X => $t[1], F => 1}); # Enqueue current SNP
	if ($#last == $opts{N}) {                   # number of SNPs in queue is N+1
	  if ($last[$#last]{X} - $last[0]{X} < $opts{W}) { # if all within window W
		map {$_->{F} = 0} @last; # all SNPs in the window of size W are "bad"
	  }
	  print STDOUT $last[0]{L} if ($last[0]{F}); # print first SNP if good
	  shift @last # dequeue first SNP
	}
  }
  # print the last few lines if applicable
  map { print $_->{L} if ($_->{F}) } @last;
}

sub pileup2fq {
  my %opts = (d=>3, D=>255, Q=>25, G=>50, l=>10);
  getopts('d:D:Q:G:l:', \%opts);
  die(qq/
Usage:   samtools.pl pileup2fq [options] <in.cns-pileup>

Options: -d INT    minimum depth        [$opts{d}]
         -D INT    maximum depth        [$opts{D}]
         -Q INT    min RMS mapQ         [$opts{Q}]
         -G INT    minimum indel score  [$opts{G}]
         -l INT    indel filter winsize [$opts{l}]
/) if (@ARGV == 0 && -t STDIN);

  my ($last_chr, $seq, $qual, @gaps, $last_pos);
  my $_Q = $opts{Q};
  my $_d = $opts{d};
  my $_D = $opts{D};

  $last_chr = '';
  while (<>) {
	my @t = split;
	if ($last_chr ne $t[0]) {
	  &p2q_post_process($last_chr, \$seq, \$qual, \@gaps, $opts{l}) if ($last_chr);
	  $last_chr = $t[0];
	  $last_pos = 0;
	  $seq = ''; $qual = '';
	  @gaps = ();
	}
	if ($t[1] - $last_pos != 1) {
	  $seq .= 'n' x ($t[1] - $last_pos - 1);
	  $qual .= '!' x ($t[1] - $last_pos - 1);
	}
	if ($t[2] eq '*') {
	  push(@gaps, $t[1]) if ($t[5] >= $opts{G});
	} else {
	  $seq .= ($t[6] >= $_Q && $t[7] >= $_d && $t[7] <= $_D)? uc($t[3]) : lc($t[3]);
	  my $q = $t[4] + 33;
	  $q = 126 if ($q > 126);
	  $qual .= chr($q);
	}
	$last_pos = $t[1];
  }
  &p2q_post_process($last_chr, \$seq, \$qual, \@gaps, $opts{l});
}

sub p2q_post_process {
  my ($chr, $seq, $qual, $gaps, $l) = @_;
  &p2q_filter_gaps($seq, $gaps, $l);
  print "\@$chr\n"; &p2q_print_str($seq);
  print "+\n"; &p2q_print_str($qual);
}

sub p2q_filter_gaps {
  my ($seq, $gaps, $l) = @_;
  for my $g (@$gaps) {
	my $x = $g > $l? $g - $l : 0;
	substr($$seq, $x, $l + $l) = lc(substr($$seq, $x, $l + $l));
  }
}

sub p2q_print_str {
  my ($s) = @_;
  my $l = length($$s);
  for (my $i = 0; $i < $l; $i += 60) {
	print substr($$s, $i, 60), "\n";
  }
}

#
# Usage
#

sub usage {
  die(qq/
Program: samtools.pl (helper script for SAMtools)
Version: $version
Contact: Heng Li <lh3\@sanger.ac.uk>\n
Usage:   samtools.pl <command> [<arguments>]\n
Command: indelFilter   filter indels generated by `pileup -c'
         snpFilter     filter SNPs generated by `pileup -c'
         pileup2fq     generate fastq from `pileup -c'
         showALEN      print alignment length (ALEN) following CIGAR
\n/);
}
