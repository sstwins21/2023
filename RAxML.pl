#!/usr/bin/perl
#running RAxML on alignment files

my ($dir, $out) = @ARGV;

@aln = `ls $dir/*.aln`;

foreach $a (@aln) {
  $a =~ s/\R//g;
  @temp = split(/\//, $a);
  print "raxmlHPC-AVX -m GTRCAT -f a -s $a -n $temp[-1] -w $out -N 2 -p 1 -x 123\n";
  #`raxmlHPC-AVX -T 2 -m GTRCAT -f d -s $a -n $temp[-1] -w $out -N 2 -p 1`;
  `raxmlHPC-AVX -m GTRCAT -f a -s $a -n $temp[-1] -w $out -N 2 -p 1 -x 123`;

}