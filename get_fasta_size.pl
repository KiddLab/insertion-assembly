#!/usr/bin/perl -w

if(@ARGV != 1) {print "one arg: fastafile\n"; exit; }
$inf=$ARGV[0];
open IN, "< $inf " or die "no in $inf\n";
$l=<IN>;
if(! ($l=~/>/))
{
  print "No '>' line?\n"; die;
}

$count=0;
while($l=<IN>)
{
  chomp($l);
  $l=~s/\s+//g;
  @bases=(split //,$l);
  $n=$#bases+1;
  $count=$count+$n;
}
close IN;
print "$count\n";