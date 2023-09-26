#!/usr/bin/perl
#counting allele freq of Clitarchus and Acanthoxyla in triplod sample


use List::MoreUtils qw(uniq);

 

#my ($vcf, $vcf2, $target_vcf) = @ARGV;
my %data;
my %ances;
my $vcf = $ARGV[0];
@cli_list = ("CLH_TP_51", "CLH_TP_52", "CLH_TP_53", "CLH_TP_54", "CLH_TP_55");
@acan_list = ("AXX_WB_10", "AXX_WB_11", "AXX_WB_13", "AXX_WB_14", "AXX_WB_9");
@inermis_1_list = ("AXP_WB_1", "AXX_WB_1", "AXX_WB_3", "AXX_WB_7");
@inermis_2_list = ("AXG_PG_11", "AXG_PG_12", "AXI_PG_1", "AXI_PG_2", "AXI_PG_3");
@prasina_1_list = ("AXP_WB_6", "AXP_WB_2021_2", "AXP_WB_2021_4", "AXP_WB_2021_8", "AXP_WB_2021_9", "R2", "R4", "R8", "R9");
@prasina_2_list = ("AXX_AB_11", "AXX_AB_13", "AXX_AB_14", "AXX_AB_15", "AXX_AB_16");



$sample_temp = `cat $vcf | grep "CHROM"`;
#$sample_temp = `grep "CHROM" "/scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/script/RNA/test/test.vcf"`;
$sample_temp =~ s/\R//g;
#print "$sample_temp\n";
@temp  = split(/\s+/, $sample_temp);

for($i = 9; $i < $#temp + 1; $i ++) {

  foreach $a (@cli_list) {
    $a =~ s/\R//g;
    #print "cli: $a\n";
    if($a eq $temp[$i]) {
      #print "cli: $a\n";
      push(@cli_sample, $i);
      last; 
    }
  }
  
   foreach $a (@acan_list) {
    $a =~ s/\R//g;
    if($a eq $temp[$i]) {
      #print "tar: $a\n";
      push(@acan_sample, $i);
      last; 
    }
  }
  
  
  foreach $a (@inermis_1_list) {
    $a =~ s/\R//g;
  
    if($a eq $temp[$i]) {
      #print "$temp[$i]\t";
      push(@inermis_1, $i);
      last; 
    }
  }
  
  foreach $a (@inermis_2_list) {
    $a =~ s/\R//g;
  
    if($a eq $temp[$i]) {
      #print "$temp[$i]\t";
      push(@inermis_2, $i);
      last; 
    }
  }
  
    foreach $a (@prasina_1_list) {
    $a =~ s/\R//g;
  
    if($a eq $temp[$i]) {
      #print "$temp[$i]\t";
      push(@prasina_1, $i);
      last; 
    }
  }
  
  foreach $a (@prasina_2_list) {
    $a =~ s/\R//g;
  
    if($a eq $temp[$i]) {
      #print "$temp[$i]\t";
      push(@prasina_2, $i);
      last; 
    }
  }
  
  
}
#print "\n";
push(@tar_sample, @inermis_1, @inermis_2, @prasina_1, @prasina_2);
foreach $a (@tar_sample) {
  print "$temp[$a]\t";
}
print "\n";

compare();





#comparing given loci with data hash to see if it present or not. printing out $pos $cli_freq $acan_freq
sub compare {

  #my ($cli_sample, $tar_sample) = @_;
  my @array = @{$_[0]};
  
  @vcf = `grep "#" -v $vcf`;
  #@vcf = `grep "#" -v "/scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/script/RNA/test/test.vcf"`;
  #print "got in: $cli_sample[0], $tar_sample[0]\n";
  foreach $a (@vcf) {
    my $PCA_loci = "";
    ($cli_freq, $tar_freq, $acan_freq, $common) = (0,0,0,0);
    $a =~ s/\R//g;
    @temp = split(/\s+/, $a);
    #print "got into $temp[0] at $temp[1]\n";
    my @alt = split(/,/, $temp[4]);
    my @cli_allele = ();
    my @acan_allele = ();
    my @tar_allele = ();
    $tar_num  = 0;
    
    
   foreach $b (@cli_sample) {
     @locus = split(/:/, $temp[$b]);
     
     if(index($locus[0], "/") != -1){
       @hap = split(/\//, $locus[0]);
       #print "$hap[0], $hap[1]\n";
     }else{
       @hap = split(/\|/, $locus[0]);
     }
     

     foreach $c (@hap) {
       if($c ne ".") {
         push(@cli_allele, $c);
         #print "$c";
       }
     }
   }
   
   foreach $b (@acan_sample) {
     @locus = split(/:/, $temp[$b]);
     
     if(index($locus[0], "/") != -1){
       @hap = split(/\//, $locus[0]);
       #print "$hap[0], $hap[1]\n";
     }else{
       @hap = split(/\|/, $locus[0]);
     }
     

     foreach $c (@hap) {
       if($c ne ".") {
         push(@acan_allele, $c);
         #print "$c";
       }
     }
   }
   
   
   
   my @unique_cli = uniq @cli_allele;
   my @unique_acan = uniq @acan_allele;
   
   #print "$#unique_cli\n";
   
   
   foreach $b (@tar_sample) {
     $sample_AD = 0;
     $needPrint = "false";
     #print "$b\n";
     @hap_AD = ();
     #($cli_freq, $tar_freq, $acan_freq, $common) = (0,0,0,0);
     @locus = split(/:/, $temp[$b]);
     if(index($locus[0], "/") != -1){
       @hap = split(/\//, $locus[0]);
       #print "$hap[0], $hap[1]\n";
     }else{
       @hap = split(/\|/, $locus[0]);
     }
     my @info = split(/\:/, $temp[8]);
    
  
     for(my $x = 0; $x < $#info + 1; $x ++) {
       if($info[$x] eq "DP") {
         $DP = $locus[$x];
       }
     }
     if($DP > 9) {
       $tar_num ++;
       for(my $x = 0; $x < $#info + 1; $x ++) {
         if($info[$x] eq "AD") {
           @temp_AD = split(/\,/, $locus[$x]);
         
           foreach $t (@temp_AD) {
             push(@hap_AD, $t/$DP);
             #print "$t/$DP\n";
           }
         
         }
       }
     }
     
     foreach $c (@hap) {
       if($c ne ".") {
         push(@tar_allele, $c);
       }else{
         $null = "null";
         $PCA_loci = "$PCA_loci$null\t";
       } 
     }
     my @unique_tar = uniq @tar_allele;
     
   
     
       
      foreach $c (@unique_tar) {
        $cli_found = "false";
        $acan_found = "false";
        $final_AD = $hap_AD[$c];
        $temp_AD = 0;
        #print "$c";
        foreach $d (@unique_cli) {
          
          if($c eq $d) {
            
            #$cli_freq = $cli_freq + $final_AD;
            #print "cli_AD: $hap_AD[$d]\n";
            $cli_found = "true";
            last
          
          }
        
        }
        
        foreach $d (@unique_acan) {
          
          if($c eq $d) {
            
            #$cli_freq = $cli_freq + $final_AD;
            #print "acan_AD: $hap_AD[$d]\n";
            $acan_found = "true";
            last
          
          }
        
        }
        
        if($acan_found eq "true" && $cli_found eq "true") {
          #print "tar_AD: $hap_AD[$d]\n";
          $common = $common + $final_AD;
        }elsif($acan_found eq "true") {
          $acan_freq = $acan_freq + $final_AD;
        }elsif($cli_found eq "true") {
          $cli_freq = $cli_freq + $final_AD;
          $sample_AD = $sample_AD + $final_AD;
        }else{
        
          $tar_freq = $tar_freq + $final_AD;
        
        }
      }
      $PCA_loci = "$PCA_loci$sample_AD\t"; 
 
       
   }
     if($common == 0 && ($#unique_acan > -1 || $tar_freq > 0 ) && $#unique_cli > -1 && $tar_num > 0) {
       $final_acan = ($acan_freq + $tar_freq);
       $final_dp = ($final_acan + $cli_freq);
       $final_acan = $final_acan/$tar_num;
       $cli_freq = $cli_freq/$tar_num;
       print "$PCA_loci\n";
       #print "$temp[0]\t$temp[1]\t$PCA_loci\n"; 
       #print "$temp[0]\t$temp[1]\t$common\t$cli_freq\t$final_acan\t$tar_num\n";
     }else{
       #print "$common\n";
     
     }
  }
}


