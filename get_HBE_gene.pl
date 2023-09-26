#!/usr/bin/perl
#adding HEB per genes
#get_HBE_gene.pl HEB.txt DAVID_ID.txt


($HEB, $DAVID) = @ARGV;

@data = `cat $HEB`;
my $results;


for $a (@data) {
  $a =~ s/\R//g;
  @temp = split(/\s+/, $a);
  $result = `grep -w $temp[0] $DAVID`;
  @temp2 = split(/\s+/, $result);
  
  if($temp2[0] ne "") {
    $temp[1] =~ s/\R//g;
    #print "$temp2[1] $temp[1]\n";
    
    if(!exists($results{$temp2[1]})){
      #print "new hash key found\n";
      $results{$temp2[1]}{hookeri} = 0;
      $results{$temp2[1]}{geisovii} = 0;
      $results{$temp2[1]}{balanced} = 0;
    }
    
    
    if($temp[1] eq "hookeri") {
      $results{$temp2[1]}{hookeri} = $results{$temp2[1]}{hookeri} + 1;
      #print "hook\n";
    }elsif($temp[1] eq "geisovii") {
      $results{$temp2[1]}{geisovii} = $results{$temp2[1]}{geisovii} + 1;
    }elsif($temp[1] eq "balanced") {
      $results{$temp2[1]}{balanced} = $results{$temp2[1]}{balanced} + 1;
    }
    
    
      
  }
}  
  

foreach my $genes (sort keys %results) {
    $total = $results{$genes}{hookeri} + $results{$genes}{geisovii} + $results{$genes}{balanced};
    
    if($results{$genes}{hookeri}/$total > 0.5) {
      print "$genes hookeri\n";
    }elsif($results{$genes}{geisovii}/$total > 0.5) {
      print "$genes geisovii\n";
    }else{
      print "$genes balances\n";
    }
}