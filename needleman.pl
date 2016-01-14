#!/usr/bin/perl
# use strict;
#use warnings;
use Switch;
# Needleman-Wunsch revisited
# input: any two input strings
# input: BLOSUM50 or matrix.txt
# output: aligned sequences 

$nl="\n";
$sp=" ";

if (@ARGV < 2) {
        print "Usage: needle.pl seq1 seq2 [options]\n";
        print "options:\n";
        print "\t m=view matrix \n";
        print "\t b=view backtrace \n";
        sys.exit();
}
$opt=();
for $i (2 .. $#ARGV) {
	$opt{$ARGV[$i]}++
}

########################## SUBROUTINES ####################################
#

# read in scoring matrix
sub getmatrix {
	my ($a) = @_;
	open(IO,$a);
	my $r=0;
	my $c=0;
	# init scoring matrix
	$sm=();
	while (my $line = <IO>) {
		chomp $line;
		if ($line =~ tr/>// or !length($line)) {
			next
		}
		@t1 = split(' ',$line);
		#print length($line),$nl;
		if (length($line) < 5) {
			$gap = $t1[0];
			close IO;
			return;
		}

		if ($r) {
			#for ($c=0;$c<@t1;$c++){
			for $c (0 .. @t1) {
				$sm{$ndx[$r-1].$ndx[$c]} = $t1[$c]
			}
		} else {
			# first row sets up the others
			for $c (0 .. @t1) {
				$ndx[$c] = $t1[$c]
			}
		}
		$r++;
	}
	close IO;
	return;
}

# read in any sequence
sub getseq {
	my ($a) = @_;
	open(IO,$a);
	my $str = "";
	while (my $line = <IO>) {
		chomp $line;
		#ignore cruft
		if ($line =~ tr/>// or !length($line)) {
			next
		}
		$protein = ($line =~ tr/ATGC//c);
		$str = $str.$line;
	}
	close IO;
	return split('',$str);
}

# return maximum by direction
sub max {
	my($a,$b,$c) = @_;
	if ($val{$a} > $val{$b}) {
		if ($val{$a} > $val{$c}) {
			return $a
		} else {
			return $c
		}
	} elsif ($val{$b} > $val{$c}) {
			return $b
		} else {
			return $c
	}
}

sub view_mat {
	if (!$opt{m}) {
		return()
	}
	print '+ ';
	for $i (0 .. $#seq1) {
		print $seq1[$i],$sp;
	}
	print $nl;
	for $i (0 .. @seq2) {
		if ($i>0) {
			print "$seq2[$i-1]:";
		}
		for $j (0 .. @seq1) {
			print $m1[$i][$j],$sp;
		}
		print $nl;
	}
}

########################## MAIN ####################################
# get all our input
#
@seq1 = getseq($ARGV[0]);
@seq2 = getseq($ARGV[1]);
if ($protein) {
	getmatrix("blosum50.txt");
} else {
	getmatrix("matrix.txt");
}
#print @seq1+0,@seq1,$nl;
#print @seq2+0,@seq2,$nl;
#print @seq1+0,$nl;
#for ($i=0;$i<@seq1;$i++) {
#	print $seq1[$i];
#}
#print $nl;
#foreach $aa (sort keys (%sm)) {
#	print "$aa $sm{$aa} $nl";
#}

# init m1 matrix gap penalties
#for ($i=0;$i<=@seq1;$i++) {
for $c (0 .. @seq1) {
	$m1[0][$c] = $c*$gap;
	$go[0][$c] = l;
}
for $r (0 .. @seq2) {
	$m1[$r][0] = $r*$gap;
	$go[$r][0] = u;
}

$val=();
# fill in calculated values
for $i (0 .. $#seq2) {
	$r = $i + 1;
	for $j (0 .. $#seq1) {
		$c = $j + 1;
		$val{u} = $m1[$r-1][$c] + $gap;
		$val{l} = $m1[$r][$c-1] + $gap;
		$val{d} = $m1[$r-1][$c-1] + $sm{$seq1[$j].$seq2[$i]};
		$dir = max(u,l,d);
		#print $val{u},$sp;
		#print $val{l},$sp;
		#print $val{d},$sp,$seq1[$j].$seq2[$i],$sp;
		#print "go $dir to $r,$c $nl";
		$m1[$r][$c] = $val{$dir};
		$go[$r][$c] = $dir;
	}
	#print $nl;
}

# see the solution matrix
view_mat;

# now the actual solution
# backtrack
#
$c = $#seq1+1;
$r = $#seq2+1;
$seq1_path = ();
$sqrt2 = ();
$seq2_path = ();

while ($c + $r > 0) {
	if ($opt{b}) {
		print "backup $r $c $go[$r][$c] $nl";
	}
	switch ($go[$r][$c]) {
	case u {
		push(@seq1_path,"-");
		push(@sqrt2," ");
		push(@seq2_path,$seq2[$r-1]);
		$r--;
	}
	case l {
		push(@seq1_path,$seq1[$c-1]);
		push(@sqrt2," ");
		push(@seq2_path,"-");
		$c--;
	}
	case d {
		push(@seq1_path,$seq1[$c-1]);
		push(@sqrt2,"|");
		push(@seq2_path,$seq2[$r-1]);
		$c--;
		$r--;
	}
	else {
		$c--;
		$r--;
	}
	}
}
print reverse(@seq1_path),$nl;
print reverse(@sqrt2),$nl;
print reverse(@seq2_path),$nl;
