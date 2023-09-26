#!/usr/bin/perl
#getting phaseblock length distribution from the results
#perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/whatshap_overlap_5.pl "/scale_wlg_nobackup/filesets/nobackup/ga02470/freebayes/lowCov/freebayes_3/whatshap/busco/all/indel_conserv/" "/scale_wlg_nobackup/filesets/nobackup/ga02470/freebayes/lowCov/freebayes_3/whatshap/busco/all/whatshap_5/" 70



($dir, $output, $filter) = @ARGV;


@folders = `ls $dir`;

######################method starts################################
$count = 0;
$geneblock_count = 0;
$unphase_port = 0;
$total_length = 0;
$gene_num = 0;
$gene_count = "false";
foreach $a (@folders) {
  $a =~ s/\R//g;
  #print "$a\n";
  @tsv = `ls $dir$a/*.tsv`;
  #print "$tsv[0]\n";
  if($#tsv > 1) {
    @consensus = ();
    #print "$a\n";
    $gene = "gene$a.fasta";
=pod
    `for files in $dir$a/*.tsv; do name="\$(cut -d "." -f 1 <<< \$files)"; grep "#" -v \$files | awk '{if(\$6 > 15) print \$2,\$4,\$5}' OFS='\t' > \$name.bed; done`;

    `bedtools intersect -a $dir$a/phased_axx.bed -b $dir$a/phased_axp.bed | bedtools intersect -a - -b $dir$a/phased_axg.bed | sort |  uniq -c | awk '{print \$2,\$3,\$4}' > $dir$a/overlap.bed`;
    `cp $dir$a/phased_axx.bed $dir$a/axx.bed`;
    `cp $dir$a/phased_axx.bed $dir$a/axp.bed`;
    `cp $dir$a/phased_axx.bed $dir$a/axg.bed`;
    `sed --regexp-extended -i s'/\:[0-9]+\-[0-9]+//' $dir$a/overlap.bed`;
    `sed --regexp-extended s'/\:[0-9]+\-[0-9]+//' $dir$a/$gene > $dir$a/renamed_gene.fasta`;
    @blocks = `perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/countN.pl $dir$a/overlap.bed $dir$a/renamed_gene.fasta`;
    
=cut  

=pod
    `for files in $dir$a/phased_*.vcf; do name=\$(echo "\$files" |rev| cut -d "_" -f 1 |rev); grep "#" \$files > $dir$a/masked_\$name; grep "#" -v \$files | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased.pl - >> $dir$a/masked_\$name; sed --regexp-extended -i s'/\:[0-9]+\-[0-9]+//' $dir$a/masked_\$name; bgzip -f $dir$a/masked_\$name; tabix -f $dir$a/masked_\$name.gz; done`;
=cut
    @blocks = `perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/countN.pl $dir$a/overlap.bed $dir$a/renamed_gene.fasta`;
    foreach $b (@blocks) {
      $b =~ s/\R//g;
      @temp = split(/\s+/, $b);
      #print "$#tsv\n";
      if($temp[3] > 300) {
        #print "$b\n";
        $count ++;
        $total_length = $total_length + $temp[3];
        $unphase_port = $unphase_port + $temp[4];
        $gene_count = "true";
        print "$temp[0]\n";
        push(@consensus, $temp[0]);

      }
    }
  
    
  }
  if($gene_count eq "true"){
    $gene_count = "false";
    $gene_num ++;
    consensus($dir, $a, $output, $filter, @consensus); 
  }
}

#$results = $total_length/$count;
#$results_2 = $unphase_port/$count;
#print "average block length that is over 300 bp is $results in $count blocks and $gene_num genes and average unphased portion is $results_2\n";




sub consensus {
  my @scaff = ();
  my ($dir, $a, $output, $filter, @scaff) = @_;
  print "going through gene $a\n";
 #`mkdir $output$a/`;
 `>$dir$a/$a.block.txt`;
  #print "############################################# $gene$i has $#scaff blocks\n";
  for($i =0; $i < $#scaff + 1; $i ++) {
    
    
    my @tempg = split(/\:/, $scaff[$i]);
    my @region = split(/\-/, $tempg[1]);
    my $fasta = `head -n1 $dir$a/gene$a.fasta`;

    my @temp2 = split(/\:/, $fasta);
    my @region2 = split(/\-/, $temp2[1]);
    my $gene_sample_name = join "", $a, "_", $i; 
    $region[0] = $region[0] + $region2[0];
    $region[1] = $region[1] + $region2[0];
    
    
    @samp_to_filter = ("$dir$a/masked_axx500.vcf.gz", "$dir$a/masked_axw2.vcf.gz", "$dir$a/masked_axp.vcf.gz", "$dir$a/masked_axg.vcf.gz", "$dir$a/masked_axx.vcf.gz");
    
    
    for($filtering = 0; $filtering < $#samp_to_filter; $filtering ++) {
      $passed = check_filter($samp_to_filter[$filtering], $dir, $a, $scaff[$i], $filter);
      if($passed eq "false") {
        last;
      }
    }
    
    if($passed eq "true") {
      `mkdir $output$a/`;
      `echo "$i $tempg[0] $region[0] $region[1]" >> $dir$a/$a.block.txt`;
      `mkdir $output$a/$i/`;
    #print "gene$a block $i has";
      dip_whatshap("$dir$a/masked_arg100.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_tp.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_pss.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_cli525.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_cli948.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_pki14.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_axw2.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_axx500.vcf.gz", $dir, $a, $scaff[$i]);
      dip_whatshap("$dir$a/masked_cli753.vcf.gz", $dir, $a, $scaff[$i]);
    
     `perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/tripvcf2cons.pl" $dir$a/renamed_gene.fasta $dir$a/masked_axx.vcf.gz $scaff[$i] "axx_$a.$i" > $output$a/$i/axx_$gene_sample_name.fasta`;
     `perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/tripvcf2cons.pl" $dir$a/renamed_gene.fasta $dir$a/masked_axg.vcf.gz $scaff[$i] "axg_$a.$i" > $output$a/$i/axg_$gene_sample_name.fasta`;
     `perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/tripvcf2cons.pl" $dir$a/renamed_gene.fasta $dir$a/masked_axp.vcf.gz $scaff[$i] "axp_$a.$i" > $output$a/$i/axp_$gene_sample_name.fasta`;


     `for files in $output$a/$i/*.fasta; do name=\$(echo "\$files" |rev| cut -d "/" -f 1 |rev | cut -d "." -f 1-4); sed -i "1s/.*/>\$name/" \$files; done`;
    }

  }


}

sub dip_whatshap {
  my ($vcf, $dir, $a, $genesblock) = @_;
  my @sample_temp = split(/\_/, $vcf);
  my @sample_name = split(/\./, $sample_temp[-1]);  
  my @tempg = split(/\:/, $genesblock);
  my @region = split(/\-/, $tempg[1]);
  $genesblock =~ s/\R//g;
  my $final_sample_name = join "", $sample_name[0], "_", $a, "_", $i; 
  
  



  my @sample_region = `cat $dir$a/filter_$sample_name[0].bed`;
  my @dip_block = ();
  #print "geneblock for $dir$a/$sample_name[0].bed: $genesblock\n";
  foreach $r (@sample_region) {
    $r =~ s/\R//g;
    @temp_region  = split(/\s+/, $r);
    #print "($temp_region[1] >= $region[0] && $temp_region[1] < $region[1]) || ($temp_region[2] > $region[0] && $temp_region[2] <= $region[1])\n";
    if(!($temp_region[1] >= $region[1]) && !($temp_region[2] <= $region[0])) {
      
      if($temp_region[1] <$region[0]) {
        $temp_region[1] = $region[0];
      }
      if($temp_region[2] > $region[1]) {
        $temp_region[2] = $region[1];
      }
      push(@dip_block, "$temp_region[1] $temp_region[2]"); 
      #print "$r: $temp_region[1] $temp_region[2]\n";
    }
    
    if($#dip_block > -1) {
      #when there are overlap
      for($c = 0; $c < $#dip_block + 1; $c ++) {
        $dip_compare = split(/\s+/, $dip_block[$c]);
        
        if($c == 0) {
          $max_dip = $dip_compare[1] - $dip_compare[0];
          $result_dip = $dip_block[$c];
          
        }elsif( $map_dip <  $dip_compare[1] - $dip_compare[0]) {
          $max_dip = $dip_compare[1] - $dip_compare[0];
          $result_dip = $dip_block[$c];
        }
      }
    } 
  }
  
  
  
  if($#dip_block > -1) {  
    #print "longest overlapping block is $result_dip\n";
    @temp_region  = split(/\s+/, $result_dip);
    
    #print "bcftools view $vcf $genesblock | grep # -v |  perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl $temp_region[0] $temp_region[1] - | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl - | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl $dir$a/renamed_gene.fasta - $genesblock $final_sample_name no > $output$a/$i/$final_sample_name.fasta\n";
    
    `bcftools view $vcf $genesblock | grep "#" -v |  perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl" $temp_region[0] $temp_region[1] - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl" - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl" $dir$a/renamed_gene.fasta - $genesblock $final_sample_name "no" > $output$a/$i/$final_sample_name.fasta`;
  }else{
    
    `bcftools view $vcf $genesblock | grep "#" -v |  perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl" -1 -1 - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl" - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl" $dir$a/renamed_gene.fasta - $genesblock $final_sample_name "no" > $output$a/$i/$final_sample_name.fasta`;
  
  }
} 


sub longranger {

  my ($vcf, $dir, $a, $genesblock) = @_;
  my $raw_fasta = "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked";
  my @blocks2 = ();
  my @sample_temp = split(/\_/, $vcf);
  my @sample_name = split(/\./, $sample_temp[-1]);

  my $final_sample_name = join "", $sample_name[0], "_", $a, "_", $i; 


  $genesblock =~ s/\R//g;
  my @tempg = split(/\:/, $genesblock);
  my @region = split(/\-/, $tempg[1]);
  my $fasta = `head -n1 $dir$a/gene$a.fasta`;

  my @temp2 = split(/\:/, $fasta);
  my @region2 = split(/\-/, $temp2[1]);

  $region[0] = $region[0] + $region2[0];
  $region[1] = $region[1] + $region2[0];

  #print "$temp[0]  $region[0] $region[1]\n";
  
 
  my $scaff1 = join "", $tempg[0], "\:", $region[0],"-", $region[1];
  my @sample = ();
  @sample = `bcftools view $vcf $scaff1 | grep "#" -v`;
  my $count = 0;
  my $blockNum = 0;
  foreach $b (@sample) {
    $b =~ s/\R//g;
    $start = "true";
    my @temp = split(/\s+/, $b);
  
  
    my @info = split(/\:/, $temp[8]);
    my $ps = 0;



      for($p = 0; $p < $#info + 1; $p ++) {
        if($info[$p] eq "PS") {
          $ps = $p;
        }
      }
  
    my @ID = split(/\:/, $temp[9]);
    if($blockNum  == 0) {
      push(@{$blocks2[$blockNum]}, $ID[$ps]); 
      push(@{$blocks2[$blockNum]}, $temp[1]);
      $blockNum ++;
        #print "before $ID[$ps]\n";
      $start = "false";
    }else{
      #print "after $blocks[-1][0] $ID[$ps] \n";
      if($blocks2[-1][0] eq $ID[$ps]) {
        push(@{$blocks2[-1]}, $temp[1]);
        #print "longer blocks $blocks[-1][0] with $ID[$ps]\n";
        #$count ++;
      }else{
        push(@{$blocks2[$blockNum]}, $ID[$ps]); 
        push(@{$blocks2[$blockNum]}, $temp[1]);
        $blockNum ++;
      }
    }

  }
  

  $count ++;
  my @longest = ();
  for($p = 0; $p < $#blocks2 + 1; $p ++) {
    if($#{$blocks2[$p]} > 1) {
      push(@longest, $p);
      #print "at index $p: $blocks2[$p][0] from $blocks2[$p][1] to $blocks2[$p][-1] which have $#{$blocks2[$p]}\n";
    }
    #print "block $i has $#{$blocks[$i]} many SNPs\n";  
  } 
  
  #when there are more than one block in overlap";
  if($#longest > 0) {
    #print "got into 1\n";
    `bcftools view $vcf $scaff1 | grep "#" -v |  perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl" $blocks2[$longest[$max]][1] $blocks2[$longest[$max]][-1] - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl" - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl" "/scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/assembly/CFSA.fasta.masked" - $scaff1 $final_sample_name "yes"> $output$a/$i/$final_sample_name.fasta`;
    
    
    
  }elsif($#longest == -1){
    #print " no AXW2 blocks\n";
    #print "got into 2\n";
    `bcftools view $vcf $scaff1 | grep "#" -v |  perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl" -1 -1 - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl" - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl" "/scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/assembly/CFSA.fasta.masked" - $scaff1 $final_sample_name "yes"> $output$a/$i/$final_sample_name.fasta`;
    
    ##when there is only one block in the overlap
  }else{
    #print "got into 3\n";
    #print "bcftools view $vcf $scaff1 | grep # -v |  perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl $blocks2[$longest[0]][1] $blocks2[$longest[0]][-1] - | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl - | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl /scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/assembly/CFSA.fasta.masked - $scaff1 $final_sample_name yes> $output$a/$i/$final_sample_name.fasta\n";
    
    `bcftools view $vcf $scaff1 | grep "#" -v |  perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl" $blocks2[$longest[0]][1] $blocks2[$longest[0]][-1] - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl" - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl" "/scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/assembly/CFSA.fasta.masked" - $scaff1 $final_sample_name "yes"> $output$a/$i/$final_sample_name.fasta`;
  
  
  }
  
 

}

sub check_filter {
  my ($vcf, $dir, $a, $genesblock, $filter) = @_;
  my @sample_temp = split(/\_/, $vcf);
  my @sample_name = split(/\./, $sample_temp[-1]);  
  my @tempg = split(/\:/, $genesblock);
  my @region = split(/\-/, $tempg[1]);
  $genesblock =~ s/\R//g;
  my $final_sample_name = join "", $sample_name[0], "_", $a, "_", $i; 
  
  


  my @sample_region = `cat $dir$a/$sample_name[0].bed`;
  my @dip_block = ();
  #print "geneblock for $dir$a/$sample_name[0].bed: $genesblock\n";
  foreach $r (@sample_region) {
    $r =~ s/\R//g;
    @temp_region  = split(/\s+/, $r);
    #print "($temp_region[1] >= $region[0] && $temp_region[1] < $region[1]) || ($temp_region[2] > $region[0] && $temp_region[2] <= $region[1])\n";
    if(!($temp_region[1] >= $region[1]) && !($temp_region[2] <= $region[0])) {
      
      if($temp_region[1] <$region[0]) {
        $temp_region[1] = $region[0];
      }
      if($temp_region[2] > $region[1]) {
        $temp_region[2] = $region[1];
      }
      push(@dip_block, "$temp_region[1] $temp_region[2]"); 
      #print "$r: $temp_region[1] $temp_region[2]\n";
    }
    
    if($#dip_block > -1) {
      #when there are overlap
      for($c = 0; $c < $#dip_block + 1; $c ++) {
        $dip_compare = split(/\s+/, $dip_block[$c]);
        
        if($c == 0) {
          $max_dip = $dip_compare[1] - $dip_compare[0];
          $result_dip = $dip_block[$c];
          
        }elsif( $map_dip <  $dip_compare[1] - $dip_compare[0]) {
          $max_dip = $dip_compare[1] - $dip_compare[0];
          $result_dip = $dip_block[$c];
        }
      }
    } 
  }
  
  
  
  if($#dip_block > -1) {  
    #print "$sample_name[0] longest overlapping block is $result_dip\n";
    @temp_region  = split(/\s+/, $result_dip);
    
    #print "bcftools view $vcf $genesblock | grep # -v |  perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl $temp_region[0] $temp_region[1] - | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl - | perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/dipvcf2cons.pl $dir$a/renamed_gene.fasta - $genesblock $final_sample_name no > $output$a/$i/$final_sample_name.fasta\n";
    
    my $block_phaserate = `bcftools view $vcf $genesblock | grep "#" -v |  perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/masking_phase.pl" $temp_region[0] $temp_region[1] - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/mask_unphased_dip.pl" - | perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/phase_rate.pl" -`;
    
    if($block_phaserate >= $filter) {
      return "true";
    }else{
      return "false";
    }
    
    #print "phaed: $sample_name[0] $block_phaserate\n";
  }else{
    #print "non_phased: $sample_name[0] 0\n";
    return "false";
  
  }
} 
