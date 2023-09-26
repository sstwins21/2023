#!/usr/bin/perl
#getting average AD in geneblock for a sample in given region

use List::MoreUtils qw(uniq);

$vcf = $ARGV[0];

$start = "false";
@samples = ();
@total = ();
@het = ();
open(my $fh, '<:encoding(UTF-8)', $vcf)
  or die "Could not open file '$vcf' $!";
 
while (my $line = <$fh>) {
  $line =~ s/\R//g;
  @temp = split(/\s+/, $line);
  
  
  if($temp[0] eq "#CHROM") {
    for($i = 9; $i < $#temp + 1; $i ++) {
    push(@samples, $temp[$i]);
    push(@AD, 0);
    push(@total, 0);
    }
   
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
      
      if($uniq[0] ne "." && $geno[$uniq[0]] ne ".") {
        $total[$i-9] = $total[$i-9] + 1;  
        if($#uniq > 0) {
          $het[$i-9] = $het[$i-9] + 1;
        }
      }
    }
  }
}

for($i = 0; $i < $#samples + 1; $i ++) {
  $rate = $het[$i]/$total[$i];
  print "$samples[$i]\t$rate\n";
}



