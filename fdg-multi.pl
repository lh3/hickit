#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (p=>0.65);
getopts('p:', \%opts);
die "Usage: fdg-multi.pl [-p probThres] <in.pairs>\n" if @ARGV == 0;

my $conf = [["2m", 0.25, 3000, 0.01, 0.2, 0.0], ["250k", 1, 1000, 0.02, 0.1, 0.0], ["50k", 1, 500, 0.03, 0.1, 5.0], ["20k", 1, 500, 0.04, 0.1, 5.0]];

my $hickit = (&dirname($0)) . '/hickit';
die 'ERROR: failed to find executable "hickit"' unless -x $hickit;

my $in = $ARGV[0];
my $prefix = $in;
$prefix =~ s/\.pairs(\.gz)?$//;

my ($prev, $next);
$next = "$prefix.$conf->[0][0].3dg.gz";
for (my $i = 0; $i < @$conf; ++$i) {
	my $cli_opts = qq/-b $conf->[$i][0] -k $conf->[$i][1] -n $conf->[$i][2] -e $conf->[$i][3] -f $conf->[$i][4] -d $conf->[$i][5]/;
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
