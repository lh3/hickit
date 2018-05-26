#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (p=>0.75, S=>1);
getopts('p:S:', \%opts);
die "Usage: fdg-multi.pl [-p probThres] [-S rngSeed] <in.phased.pairs>\n" if @ARGV == 0;

my $conf = [["4m",   0.2, 20000, 0.02, 0.0],
			["1m",   0.2,  5000, 0.02, 0.0],
			["200k", 0.2,  1500, 0.02, 0.0],
			["50k",  0.2,   750, 0.02, 5.0],
			["20k",  0.2,   750, 0.02, 5.0]];

my $hickit = (&dirname($0)) . '/hickit';
die 'ERROR: failed to find executable "hickit"' unless -x $hickit;

my $in = $ARGV[0];
my $prefix = $in;
$prefix =~ s/\.pairs(\.gz)?$//;

my ($prev, $next);
for (my $i = 0; $i < @$conf; ++$i) {
	my $cli_opts = qq/-S $opts{S} -b $conf->[$i][0] -k $conf->[$i][1] -n $conf->[$i][2] -e $conf->[$i][3] -d $conf->[$i][4]/;
	$next = $i == @$conf - 1? "$prefix.3dg.gz" : "$prefix.$conf->[$i][0].3dg.gz";
	if ($i == 0) {
		print "$hickit bin -g -p $opts{p} $cli_opts $in 2> $prefix.3dg.log | gzip > $next\n";
	} else {
		print "$hickit bin -g -p $opts{p} -i $prev $cli_opts $in 2>> $prefix.3dg.log | gzip > $next\n";
	}
	$prev = $next;
}

sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}

sub make_cli_opts {
	my $c = shift(@_);
	my $i = shift(@_);
	return qq/-b $conf->[$i][0] -k $conf->[$i][1] -n $conf->[$i][2] -e $conf->[$i][3]/;
}
