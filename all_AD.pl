#!/usr/bin/perl
#getting average AD in geneblock for a sample in given region

use List::MoreUtils qw(uniq);

$vcf = $ARGV[0];

$start = "false";
@samples = ();

open(my $fh, '<:encoding(UTF-8)', $vcf)
  or die "Could not open file '$vcf' $!";
 
while (my $line = <$fh>) {
  $line =~ s/\R//g;
  @temp = split(/\s+/, $line);
  
  
  if($temp[0] eq "#CHROM") {
    for($i = 9; $i < $#temp + 1; $i ++) {
      print "$temp[$i]\t";
    }
    print "\n";
    $start = "true";  
    
  }elsif($start eq "true") {  
    @geno = ();
    
    push(@geno, $temp[3]);
    @alt = split(/\,/, $temp[4]);
    foreach $g (@alt) {
      push(@geno, $g);
    } 
    
    for($i = 9; $i < $#temp + 1; $i ++) {
      @info = split(/\:/, $temp[$i]);
      
      if(index($info[0], "/") != -1) {
        @hap = split(/\//, $info[0]); 
      }else{
        @hap = split(/\|/, $info[0]);
      }
      
      @uniq = uniq(@hap);
      
      if($uniq[0] eq "." || $geno[$uniq[0]] eq ".") {
        print "NA\t";
      }else{
        print "$info[2]\t";
      }
    }
    
    print "\n";
  
  }
}



