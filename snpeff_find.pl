#!/usr/bin/perl
#getting SNPs from gene list and snpeff annotation
use List::MoreUtils qw(uniq);


($species, $name) = @ARGV;


@region = ();
%annot = ();
#getting sample index from header

@cli_list = ("CLH_TP_51", "CLH_TP_52", "CLH_TP_53", "CLH_TP_54", "CLH_TP_55");
@acan_list = ("AXX_WB_10", "AXX_WB_11", "AXX_WB_13", "AXX_WB_14", "AXX_WB_9");
  
if($species eq "cli") {
  $path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_cli/new_2/cor/sorted_", $name, ".bed";
  $vcf = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/cli.ann.vcf";
  #$vcf = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter4/test/tset.vcf";
  $tar_bed = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/cli.bed";

}elsif($species eq "acan") {
  $path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_acan/new_2/cor/sorted_", $name, ".bed";;
  $vcf = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/axw2.ann.vcf";
  $tar_bed = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/acan.bed";

}elsif($species eq "arg100") {
  $path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/arg100/HEB_2/cor/sorted_", $name, ".bed";;
  $vcf = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/arg100.ann.vcf";
  $tar_bed = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/arg100.bed";

}

@gene_list = `cat /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/$species/$name.david`;  
  
if($name eq "inermis_1") {
  @sample_list = ("AXP_WB_1", "AXX_WB_1", "AXX_WB_3", "AXX_WB_7");
}elsif($name eq "inermis_2") {
  @sample_list = ("AXG_PG_11", "AXG_PG_12", "AXI_PG_1", "AXI_PG_2", "AXI_PG_3");
}elsif($name eq "prasina_1") {
  @sample_list = ("AXP_WB_6", "AXP_WB_2021_2", "AXP_WB_2021_4", "AXP_WB_2021_8", "AXP_WB_2021_9", "R2", "R4", "R8", "R9");
}elsif($name eq "prasina_2") {
  @sample_list = ("AXX_AB_11", "AXX_AB_13", "AXX_AB_14", "AXX_AB_15", "AXX_AB_16");
}else{
  #print "wrong sample name \n";
}


$sample_temp = `grep "CHROM" "$vcf"`;
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
      #print "acan: $a\n";
      push(@acan_sample, $i);
      last; 
    }
  }
  
  
  foreach $a (@sample_list) {
    $a =~ s/\R//g;
    if($a eq $temp[$i]) {
      #print "tar: $a\n";
      push(@tar_sample, $i);
      last; 
    }
  }
}


#getting annotations from snpeff vcf file
@vcf = `grep "#" -v "$vcf"`;

foreach $a (@vcf) {
  $a =~ s/\R//g;
  @temp = split(/\s+/, $a);
  
  $annot{$temp[0]}{$temp[1]} = $a;
  
}


foreach $a (@gene_list) {
   $a =~ s/\R//g;
   @gene_temp = split(/\s+/, $a);
   $bed = `grep -w @gene_temp[0] $tar_bed`;
   $bed =~ s/\R//g;
   
   if($bed eq "") {
     print "somthing is wrong at getting gene region\n";
   }else{
     #print "$bed\n";
     @temp = split(/\s+/, $bed);
     @tar_region = `awk '{if("$temp[0]" == \$1 && \$2 >= $temp[1] && \$2 <= $temp[2]) print \$1, \$2, \$4}' $path`;
     
     foreach $b (@tar_region) {
       $b =~ s/\R//g;
       @annot_tar = ();
       @temp2 = split(/\s+/, $b);
       if(exists $annot{$temp2[0]}{$temp2[1]} ) {
         @temp3 = split(/\s+/,$annot{$temp2[0]}{$temp2[1]});
         print "$temp2[0]\t$temp2[1]\t";
         my @alt = split(/,/, $temp3[4]);
         my @cli_allele = ();
         my @acan_allele = ();
         my @tar_allele = ();
         
         my @annot_temp = split(/;/, $temp3[7]);
         for($i = 0; $i < $#annot_temp + 1; $i ++) {
           if(substr($annot_temp[$i], 0, 3) eq "ANN") {
             @annot = split(/,/, $annot_temp[$i]);
             last;
             
           }
         }     
         #print "annotation: $annot_temp[-1]\n";
         #my @annot = split(/,/, $annot_temp[-1]);
       
         push(@annot_tar, "ref");
         foreach $c (@annot) {
           @temp4 = split(/\|/, $c);
           #print "test: $c \n";
           push(@annot_tar, $temp4[1]);
         }
         #print "\n";
    
    
         foreach $c (@cli_sample) {
           @locus = split(/:/, $temp3[$c]);
     
           if(index($locus[0], "/") != -1){
             @hap = split(/\//, $locus[0]);
           #print "$hap[0], $hap[1]\n";
           }else{
             @hap = split(/\|/, $locus[0]);
           }
     
        
           foreach $d (@hap) {
             if($d ne ".") {
               push(@cli_allele, $d);
             }
           }
        }
   
        foreach $c (@acan_sample) {
          #print "acan: $c\n";
          @locus = split(/:/, $temp3[$c]);
      
          if(index($locus[0], "/") != -1){
            @hap = split(/\//, $locus[0]);
          #print "$hap[0], $hap[1]\n";
          }else{
            @hap = split(/\|/, $locus[0]);
          }  
      
 
          foreach $d (@hap) {
       
            if($d ne ".") {
              push(@acan_allele, $d);
            }
          }
        }
    
    
    
        my @unique_cli = uniq @cli_allele;
        my @unique_acan = uniq @acan_allele;
        #print " cli: ";
        foreach $c (@unique_cli) {
          #print " $c";
        }
        #print " acan: ";
        foreach $c (@unique_acan) {
          #print " $c";
        }
      
        foreach $c (@tar_sample) {
          @locus = split(/:/, $temp3[$c]);
     
         if(index($locus[0], "/") != -1){
           @hap = split(/\//, $locus[0]);
           #print "$hap[0], $hap[1]\n";
         }else{
           @hap = split(/\|/, $locus[0]);
         }
         foreach $c (@hap) {
           if($c ne ".") {
             push(@tar_allele, $c);
           } 
         }
        }
        my @unique_tar = uniq @tar_allele;
      
        @cli_annot = ();
        @acan_annot = ();

        for($u = 0; $u < $#unique_tar + 1; $u ++) {
          if($u == 0) {
            print "$annot_tar[$u]";
          }else{
            print ",$annot_tar[$u]";
          }
        }
        print "\t$temp2[2]\t$a\n";     
      }
    }

     
   }

}
