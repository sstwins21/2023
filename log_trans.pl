#!/usr/bin/perl
#transform the output file from snpeff_find.pl to R database
use List::MoreUtils qw(uniq);

($species, $name) = @ARGV;

@temp = `awk '{print \$5}' "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/$species/$name.snpeff"`;
@genes = uniq @temp;
@rates = ();
@HEB = ();
@total = ();
@FPKM = ();

if($species eq "cli") {
  $path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_cli/new_2/cor/", $name, ".HBE2.txt";
  $david = "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/david_2.txt";
  $fpkm_path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_cli/new_2/cor/", $name, ".HEB";
}elsif($species eq "acan") {
  $path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_acan/new_2/cor/", $name, ".HBE2.txt";
  $david = "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/acan_david_2.txt";
  $fpkm_path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_acan/new_2/cor/", $name, ".HEB";
}elsif($species eq "arg100") {
  $path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/arg100/HEB_2/cor/", $name, ".HBE2.txt";
  $david = "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/arg100_david_2.txt";
  $fpkm_path = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/arg100/HEB_2/cor/", $name, ".HEB";
}

#going trhough gene list
foreach $a (@genes) {
  chomp $a;
  ($syn_num, $non_num) = (0,0);
  
  #print "$a\n";
  @current = `grep -w "$a" "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/annotation/$species/$name.snpeff"`;
  #chomp $current[0];
  
  #looking at each gene
  foreach $b (@current) {
    chomp $b;
    #print "$b\n";
    $syn = "false";
    $non = "false";
    $dir = "";
    
    @snpeff_temp  = split(/\s+/, $b);
    if($snpeff_temp[3] eq "hookeri") {
      $dir = "cli";
    }elsif($snpeff_temp[3] eq "geisovii") {
      $dir = "acan";
    }
    @temp2  = split(/\,/, $snpeff_temp[2]);
      #looping through the annotation of the locus
    foreach $c (@temp2) {
      chomp $c;
      if($c ne "intergenic_region") {
        if($c eq "synonymous_variant" || $c eq "start_retained_variant" || $c eq "stop_retained_variant" ) {
          $syn = "true";
        }elsif($c eq "missense_variant") {
          $non = "true";
        }
      }
    }
      
    if($non eq "true") {
      $non_num = $non_num + 1;
    }elsif($syn eq "true") {
      $syn_num = $syn_num + 1;
    }
    
  }
  
  if($syn_num > 0 && $non_num > 0 && "balanced" ne $snpeff_temp[3] && $dir eq $snpeff_temp[-1]) {
  #if($syn_num > 0 && $non_num > 0 && $dir eq $snpeff_temp[-1]) {
  #if($syn_num > 0 && $non_num > 0) {
    $rate = $non_num/$syn_num;
    $HBE_temp = `grep -w "$a" $path`;
    @temp  = split(/\s+/, $HBE_temp);
    if($temp[5] ne "") {
      push(@rates, $rate);
      chomp $temp[5];
      push(@HEB, $temp[5]);
      push(@total, $syn_num + $non_num);
      
      $ID = `grep -w "$a" $david | awk '{print \$2}'`;
      chomp $ID;
      #print "$a $ID\n";
      $fpkm = `grep -w "$ID" $fpkm_path`;
      chomp $fpkm; 
      #print "$a $ID $fpkm\n";
      @fpkm_temp = split(/\s+/, $fpkm);
      if($fpkm_temp[-1] eq "geisovii") {
        
        push(@FPKM, $fpkm_temp[2]);
      }elsif($fpkm_temp[-1] eq "hookeri") {
       
        push(@FPKM, $fpkm_temp[1]);
      }else{
        print "wrong at FPKM step\n";
      }
    }
  }
  #$HEB = `grep -w "$a" 
  #push 
}

for($i = 0; $i < $#rates + 1; $i ++) {
  print "$FPKM[$i]\t$HEB[$i]\t$rates[$i]\t$total[$i]\n";
  
  
}