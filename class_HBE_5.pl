#!/usr/bin/perl
#classifing if the gene is biased/ mono / balanced
#perl /scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/script/RNA/class_HBE.pl

$bed_file = $ARGV[0];


if($ARGV[1] eq "cli") {
  @gene_bed = `cat "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/cli_gene.bed"`;
}elsif($ARGV[1] eq "acan") {
  @gene_bed = `cat "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/acan_gene.bed"`;
}elsif($ARGV[1] eq "arg100") {
  @gene_bed = `cat "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/arg100/arg100_gene.bed"`;
}
foreach $a (@gene_bed) {
  $a =~ s/\R//g;
  @temp = split(/\s+/, $a);
  #@name = split(/=/, $temp[3]);
  #print "$temp[3] $name[1]\n";
  
  @data = `awk '{if(\$1 == "$temp[0]" && \$2 >= $temp[1] && \$3 <= $temp[2]) print \$0}' $bed_file`;
  #print "$name[1] $temp[0] $temp[1] $temp[2]: $#data\n";

  if($#data > -1) {
    #print "$name[1] $temp[0] $temp[1] $temp[2]: $#data\n";
    ($final_acan, $final_cli, $final_BHE, $loci_acan, $loci_cli, $loci_BHE) = (0,0,0,0,0,0);
    $prev = "";
    @cli_dos = ();
    @acan_dos = ();
    
    
    
    foreach $b (@data) {
      $b =~ s/\R//g;
      @temp2 = split(/\s+/, $b);
      #print "$b\n";
      #print "$prev ne $temp2[2] && ($loci_cli > 0 || $loci_BHE > 0 || $loci_acan > 0\n";
      if($prev ne $temp2[2] && ($loci_cli > 0 || $loci_BHE > 0 || $loci_acan > 0)) {
        
        $total = $loci_cli + $loci_acan + $loci_BHE;
        #print "$prev ne $temp2[2] && ($loci_cli > 0 || $loci_BHE > 0 || $loci_acan > 0\n";
        if($loci_cli/$total > 0.7) {
          $final_cli ++;
        }elsif($loci_acan/$total > 0.7) {
          $final_acan ++;
        }else{
          $final_BHE ++;
          
        }
        ($loci_acan, $loci_cli, $loci_BHE) = (0,0,0);
      }
      
      if($temp2[3] eq "geisovii") {
        $loci_acan ++;
        push(@acan_dos, $temp2[4]);
      }elsif($temp2[3] eq "hookeri") {
        $loci_cli ++;
        push(@cli_dos, $temp2[4]);
      }elsif($temp2[3] eq "balanced") {
        #print "got in to adding\n";
        $loci_BHE ++; 
      }
      $prev = $temp2[2];
      #print "$prev\n";
    }
    
    #$total = $loci_cli + $loci_acan + $loci_BHE;
    $total = $loci_cli + $loci_acan;
    if($total == 0) {
      $total = 1;
    }
    #print "$prev ne $temp2[2] && ($loci_cli > 0 || $loci_BHE > 0 || $loci_acan > 0\n";
    if($loci_cli/$total > 0.7) {
      $final_cli ++;
    }elsif($loci_acan/$total > 0.7) {
      $final_acan ++;
    }else{
      $final_BHE ++;
      
    }

    #$total = $final_cli + $final_acan + $final_BHE;
    $total = $final_cli + $final_acan;
    if($total == 0) {
      $total = 1;
    }
    $my_sum = 0;
    #print "$name[1] $final_cli $final_acan $final_BHE $#data\n";
    #print "$name[1]\t";
    print "$temp[3]\t";
    if($final_cli/$total > 0.7) {
      print "hookeri\t";
      foreach $a (@cli_dos) {
        $my_sum = $my_sum + $a;
      }
      $average = $my_sum/($#cli_dos + 1);
    }elsif($final_acan/$total > 0.7) {
      print "geisovii\t";
      foreach $a (@acan_dos) {
        $my_sum = $my_sum + $a;
      }
      $average = $my_sum/($#acan_dos + 1);
    }else{
      print "balanced\t";
      $average = 0;  
    }
    
    
    #print "$average\n";
    print "$final_cli\t$final_acan\t$final_BHE\t$average\n";

  }
}

