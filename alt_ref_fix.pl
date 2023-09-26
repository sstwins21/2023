#!/usr/bin/perl
#Removing alt ref to account more SNPs from whatshap
 

use List::MoreUtils qw(uniq);



open(my $in, '<:encoding(UTF-8)', $ARGV[0])
or die "Could not open file '$vcf' $!";
 
while (my $line = <$in>) {
  chomp $line;
  
  if(substr($line, 0, 1,) eq "#") {
    print "$line\n";
  
  }else{
  
    @temp = split(/\s+/, $line);

   
    @info = split(/\:/, $temp[9]);
    if(index($info[0], "/") != -1) {
      @hap = split(/\//, $info[0]); 
      $sep="\/";
    }else{
      @hap = split(/\|/, $info[0]);
      $sep="\|";
    }
    #print "$sep\n";
    @alt = split(/\,/, $temp[4]);
    @ad = split(/\,/, $info[1]);
    #my @hap =  split(/\//, $info[0]);
    my @uniq = uniq @hap;
    #print "@uniq ";
    #print "@ad\n";
    #print "@alt $#alt ";

    if($uniq[-1] > 1 ) {
      #print "##########run ";
      $start_geno = 1;
      %hash_geno = ();
      %hash_alt = ();
      @new_geno = ();
      @new_ad = ();
      @new_alt = ();
      push(@new_ad, 0);
      foreach $a (@hap) {
        if($a == 0) {
          push(@new_geno, $a);
          
        }else{
          
          if(!exists($hash_geno{$a})){
            $hash_geno{$a} = $start_geno;
            $start_geno = $start_geno + 1;
            push(@new_alt, $alt[$a-1]);
            push(@new_ad, $a);
          }
          push(@new_geno, $hash_geno{$a});
          
        }
      
      }
      #print "@new_geno @new_alt @new_ad";
      for($i = 0; $i < $#temp + 1; $i ++) {
        if($i == 4) {
          for($o = 0; $o < $#new_alt + 1; $o ++) {
            print "$new_alt[$o]";
            if($o == $#new_alt) {
              print "\t";
            }else{
              print ",";
            }
          }
        
        }elsif($i == 9) {
          for($o = 0; $o < $#new_geno + 1; $o ++) {
            print "$new_geno[$o]";
            if($o == $#new_geno) {
              print ":";
            }else{
              print "$sep";
            }
          }
          
          for($o = 0; $o < $#new_ad + 1; $o ++) {
            print "$ad[$new_ad[$o]]";
            if($o == $#new_ad) {
              print ":";
            }else{
              print ",";
            }
          }
          
          for($o = 2; $o < $#info + 1; $o ++) {
            print "$info[$o]";
            if($o == $#info) {
              print "\n";
            }else{
              print ":";
            }
          }
        
        }else{
          print "$temp[$i]\t";
        }
      
      
      }
      
    }else{
      print "$line\n";
    }
    
    
    #print "\n";
    
    
    
  }
}