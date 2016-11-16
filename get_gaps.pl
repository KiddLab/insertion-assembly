#!/usr/bin/perl -w

if(@ARGV != 1) {print "1 arg, fasta file, will return positions of Ns\n"; exit; }

$inf=$ARGV[0];

$cmd="get_fasta_size.pl $inf";
$size=`$cmd`;
chomp($size);

@arry=();

$str="";
open IN, "< $inf " or die "NO IN";
$l=<IN>;
chomp($l);
$name=$l;
$name=~s/>//;
$l=<IN>;
chomp($l);
$l=~s/\s+//g;
$str=$l;
while($l=<IN>)
{
 chomp($l);
 $l=~s/\s+//g;
 $str.=$l;
}
close IN;

@arry=(split //,$str);
unshift(@arry, "?");


$lp=0;
$rp=0;
#step through getting regions.
$i=1;
$tg=0;
while($i <=$size)
{
 while($i <=$size and ( uc($arry[$i]) ne "N") )  {$i++; }
 #in reg start
 $lp=$i;
 while($i <= $size and ( uc($arry[$i]) eq "N") ) {$i++; }
 $rp=$i-1; #got rp
 if($lp <= $size)
 {
  print "$name\t$lp\t$rp\tGAP\n"; 
  $tg=$tg+($rp-$lp)+1;
 }
} #end while i < re
#print "Total gaps=$tg\n";