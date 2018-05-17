#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (k=>1, r=>2.0, p=>0.65);
getopts('k:r:p:', \%opts);
die "Usage: fdg-multi.pl [-k repCoef] [-r refRadius] [-p probThres] <in.pairs>\n" if @ARGV == 0;

my $conf = [["2m", 2000, 0.01], ["250k", 1000, 0.02], ["50k", 500, 0.03], ["20k", 500, 0.04]];

my $hickit = (&dirname($0)) . '/hickit';
die 'ERROR: failed to find executable "hickit"' unless -x $hickit;

my $in = $ARGV[0];
my $prefix = $in;
$prefix =~ s/\.pairs(\.gz)?$//;

my ($prev, $next);
$next = "$prefix.$conf->[0][0].3dg.gz";
print "$hickit bin -g -p $opts{p} -k $opts{k} -r $opts{r} -b $conf->[0][0] -n $conf->[0][1] -e $conf->[0][2] $in 2> $prefix.3dg.log | gzip > $next\n";
$prev = $next;
for (my $i = 1; $i < @$conf; ++$i) {
	$next = $i == @$conf - 1? "$prefix.3dg.gz" : "$prefix.$conf->[$i][0].3dg.gz";
	print "$hickit bin -g -p $opts{p} -k $opts{k} -r $opts{r} -i $prev -b $conf->[$i][0] -n $conf->[$i][1] -e $conf->[$i][2] $in 2>> $prefix.3dg.log | gzip > $next\n";
	$prev = $next;
}

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
