#!/usr/bin/perl -w
=head1 NAME        

createSegFilesFromAPOLLOH.pl

-head1 SYNOPSIS

=head1 OPTIONS 
    
    -calls                use if want calls in results
    -outfile|o <string>
    -infile|i <string>
  
=head1 DESCRIPTION

=head1 CONTACT

Gavin Ha <gha@bccrc.ca>

=cut

use strict;
use List::Util qw[min max];
use File::Basename;
use Getopt::Long;

sub usage () {
    exec('perldoc', $0);
    exit;
}

my ($outfile, $infile, $calls, $help);

GetOptions (
	    'infile|i=s' => \$infile,	    
	    'calls' => \$calls,
            'outfile|o=s' => \$outfile,
            'help|?' => \$help
            );

if($help) {
    &usage();
}
print "Parameters:\ninfile=$infile\noutfile=$outfile\ncalls=$calls\n";

my ($name,$path,$suffix) = fileparse($infile);
my ($id,@jnk) = split(/\_/,$name);
my $nameOut = $id;

open OUTFILE, ">$outfile\.tmp";
print OUTFILE "Sample\tChromosome\tStart Position(bp)\tEnd Position(bp)\tLength(Mb)\tMedian_Ratio\tSeg_CN\n";
open(SEGFILE, $infile) || die("Can't open $infile!\n");
my $line = <SEGFILE>; chomp($line);
my($chr, $start, $ref, $nRef, $N, $ratio, $cn, $state, $call) = split(/\t/,$line);
my $end = $start;
my @totalRatio=();
push(@totalRatio,max($ratio,1-$ratio));
while($line=<SEGFILE>) {
chomp($line);
my ($chrS, $startS, $refS, $nRefS, $NS, $ratioS, $cnS, $stateS, $callS) = split(/\t/,$line);	
if ($chrS!=$chr || $callS ne $call){
	my $medianRatio = max($ratio,1-$ratio); #1 position only
	$medianRatio = median(@totalRatio) if (($end-$start+1)>1); #more than 1 position
	my $output = $nameOut . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . ($end-$start+1). "\t" . sprintf("%.4f",$medianRatio);	
	$output = $output . "\t" .  $call if ($calls);
	print OUTFILE $output . "\n";
	#reset
	@totalRatio = ();
	$start = $startS; 
}else{
	
	#push(@totalRatio,max($ratioS,1-$ratioS));
	#assign current state to previous variables
	#new end, but start still the same
}
($chr, $end, $ref, $nRef, $N, $ratio, $cn, $state, $call) = ($chrS, $startS, $refS, $nRefS, $NS, $ratioS, $cnS, $stateS, $callS);
push(@totalRatio,max($ratioS,1-$ratioS));
}

#final output
my $medianRatio = max($ratio,1-$ratio); #1 position only
$medianRatio = median(@totalRatio) if (($end-$start+1)>1); #more than 1 position
my $output = $nameOut . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . ($end-$start+1). "\t" . sprintf("%.4f",$medianRatio);	
$output = $output . "\t" .  $call if ($calls);
print OUTFILE $output . "\n";
close OUTFILE;
my $returnCode = `sort -k2,2n -k3,3n $outfile\.tmp > $outfile; rm $outfile\.tmp;`;
close SEGFILE;


sub median { 
	my (@array_ref) = @_; 
	my $count = scalar @array_ref; 
	# Sort a COPY of the array, leaving the original untouched 
	my @array = sort { $a <=> $b } @array_ref; 
	if ($count % 2) { 
		return $array[int($count/2)]; 
	} else { 
		return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
} 

#Read more: http://wiki.answers.com/Q/How_can_you_calculate_the_average_and_median_in_perl_by_subroutine#ixzz1LcMWdXEE
#sub median{
#	my @a = sort @_;
#	return ($a[$#a/2] + $a[@a/2]) / 2;
#}
