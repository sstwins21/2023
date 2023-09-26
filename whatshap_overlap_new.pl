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




    print "compress begins\n";
    `for files in $dir$a/phased_*.vcf; do name=\$(echo "\$files" |rev| cut -d "_" -f 1 |rev);  echo "\$name" ;bgzip -c \$files > $dir$a/\$name.gz; tabix -f $dir$a/\$name.gz; done`;
    print "compress done\n";






    @blocks = `perl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/countN.pl $dir$a/overlap.bed`;
    
    foreach $b (@blocks) {
      $b =~ s/\R//g;
      @temp = split(/\s+/, $b);
      #print "$#tsv\n";
      #print "temp $temp[3]\n";
      if($temp[3] > 300) {
        #print "$b\n";
        $count ++;
        $total_length = $total_length + $temp[3];
        $unphase_port = $unphase_port + $temp[4];
        $gene_count = "true";
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
    
    
    @samp_to_filter = ("$dir$a/axx500.vcf.gz", "$dir$a/axw2.vcf.gz", "$dir$a/axp.vcf.gz", "$dir$a/axg.vcf.gz", "$dir$a/axx.vcf.gz");
    @all_samp_to_filter = ("$dir$a/arg100.vcf.gz", "$dir$a/axg.vcf.gz", "$dir$a/axp.vcf.gz", "$dir$a/axw2.vcf.gz", "$dir$a/axx500.vcf.gz", "$dir$a/axx.vcf.gz", "$dir$a/cli525.vcf.gz", "$dir$a/cli753.vcf.gz", "$dir$a/cli948.vcf.gz", "$dir$a/pki14.vcf.gz", "$dir$a/pss.vcf.gz", "$dir$a/tp.vcf.gz","$dir$a/AXW2.vcf.gz");
    
    for($filtering = 0; $filtering < $#samp_to_filter + 1; $filtering ++) {
      $passed = check_filter($samp_to_filter[$filtering], $dir, $a, $scaff[$i], $filter);
      #print "is it passed $passed\n";
      if($passed eq "false") {
        last;
      }
    }
    
    if($passed eq "true") {
      for($filtering = 0; $filtering < $#all_samp_to_filter + 1 ; $filtering ++) {
        $passed = check_filter_2($all_samp_to_filter[$filtering], $dir, $a, $scaff[$i], $filter);
    
        if($passed eq "false") {
          last;
        }
      }
    }  
    
    if($passed eq "true") {
      `mkdir $output$a/`;
      `echo "$i $scaff[$i]" >> $dir$a/$a.block.txt`;
      `mkdir $output$a/$i/`;
      #print "gene$a block $i has";
      #print "$scaff[$i]\n";
      
      phasing("$dir$a/arg100.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/tp.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/pss.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/cli525.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/cli948.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/pki14.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/axw2.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/axx500.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/cli753.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/AXW2.vcf.gz", $dir, $a, $scaff[$i], "dip");
      phasing("$dir$a/axx.vcf.gz", $dir, $a, $scaff[$i], "trip");
      phasing("$dir$a/axp.vcf.gz", $dir, $a, $scaff[$i], "trip");
      phasing("$dir$a/axg.vcf.gz", $dir, $a, $scaff[$i], "trip");
      
      
      
      
      
     #`for files in $output$a/$i/*.fasta; do name=\$(echo "\$files" |rev| cut -d "/" -f 1 |rev | cut -d "." -f 1-4); sed -i "1s/.*/>\$name/" \$files; done`;
    }

  }


}


sub check_filter {
  my ($vcf, $dir, $a, $genesblock, $filter) = @_;
  my @sample_temp = split(/\//, $vcf);
  my @sample_name = split(/\./, $sample_temp[-1]);  
  my @tempg = split(/\:/, $genesblock);
  my @region = split(/\-/, $tempg[1]);
  $genesblock =~ s/\R//g;
  my $final_sample_name = join "", $sample_name[0], "_", $a, "_", $i; 
  
  #print "$sample_name[0]\n";


  my @sample_region = `cat $dir$a/phased_$sample_name[0].bed`;
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

sub phasing {
  my ($vcf, $dir, $a, $genesblock, $ploidy) = @_;
  my @sample_temp = split(/\//, $vcf);
  my @sample_name = split(/\./, $sample_temp[-1]);  
  my @tempg = split(/\:/, $genesblock);
  my @region = split(/\-/, $tempg[1]);
  $genesblock =~ s/\R//g;
  my $final_sample_name = join "", $sample_name[0], "_", $a, "_", $i; 
  
  
  print "$geneblock\n";

  my @sample_region = `bcftools view -H -r $genesblock $vcf`;
  my %phase_ID = ();
  my $max = 0;
  $max_ID = "failed";
  foreach $p (@sample_region) {
    $p =~ s/\R//g;
    @Ptemp = split(/\s+/, $p);
    my @p_ID = split(/:/, $Ptemp[9]);
  
    if(index($p_ID[0], "/") != -1){
      @hap = split(/\//, $p_ID[0]);
      #print "$hap[0], $hap[1]\n";
    }else{
      @hap = split(/\|/, $p_ID[0]);
    }    
    
    
    
    
    if($hap[0] ne "." && $p_ID[-1] ne "." && index($p_ID[-1], ",") == -1) {
      #print "got in \n";
      if(exists($phase_ID{$p_ID[-1]})) {
        $phase_ID{$p_ID[-1]} = $phase_ID{$p_ID[-1]} + 1
        
      }else{
        $phase_ID{$p_ID[-1]} = 1;
      }
      
      if($max < $phase_ID{$p_ID[-1]}) {
        $max = $phase_ID{$p_ID[-1]};
        $max_ID = "$p_ID[-1]";
      }
      
    }
    
  }
  
  
  `bcftools view -h -r $genesblock $vcf > $dir$a/$final_sample_name.vcf`;
  unless(open FILE, '>>', "$dir$a/$final_sample_name.vcf") {
	  # Die with error message 
	  # if we can't open it.
  	die "nUnable to create $name/n";
  }
  foreach $p (@sample_region) {
    $p =~ s/\R//g;
    @Ptemp = split(/\s+/, $p);
    my @p_ID = split(/:/, $Ptemp[9]);
    if(index($p_ID[0], "/") != -1){
      @hap = split(/\//, $p_ID[0]);
      #print "$hap[0], $hap[1]\n";
    }else{
      @hap = split(/\|/, $p_ID[0]);
    }      
    
    
    
    if($p_ID[-1] eq $max_ID) {
      print FILE "$p";  
    
    }else{
      for($q = 0; $q < $#Ptemp; $q ++) {
        print FILE "$Ptemp[$q]\t";
      
      }
      if($ploidy eq "dip") { 
        print FILE "$hap[0]/$hap[1]";
      }elsif($ploidy eq "trip") {
        print FILE "$hap[0]/$hap[1]/$hap[2]";
      }
      for($q = 0; $q < $#p_ID; $q ++) {
        print FILE ":$q";
      }
    }
    print FILE "\n";
  
  
  
  
  }
  close(FILE);
  
  
  `bgzip -f $dir$a/$final_sample_name.vcf`;
  `tabix -f $dir$a/$final_sample_name.vcf.gz`;
  
  if($ploidy eq "dip") {
    `samtools faidx /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked $genesblock | bcftools consensus -H 1pIu $dir$a/$final_sample_name.vcf.gz > $output$a/$i/$final_sample_name.fasta`;
    `sed -i 's/$genesblock/$final_sample_name.1/g' $output$a/$i/$final_sample_name.fasta`;
    `samtools faidx /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked $genesblock | bcftools consensus -H 2pIu $dir$a/$final_sample_name.vcf.gz >> $output$a/$i/$final_sample_name.fasta`;
    `sed -i 's/$genesblock/$final_sample_name.2/g' $output$a/$i/$final_sample_name.fasta`;
    
  }elsif($ploidy eq "trip") {
    
    `perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/PhD/Chapter3/tripvcf2cons.pl" /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/whatshap/scaffolds/$tempg[0].fasta $dir$a/$final_sample_name.vcf.gz $genesblock $final_sample_name > $output$a/$i/$final_sample_name.fasta`;
    #print "perl tripvcf2cons.pl /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/whatshap/scaffolds/$tempg[0].fasta $dir$a/$final_sample_name.vcf.gz $genesblock $final_sample_name > $output$a/$i/$final_sample_name.fasta\n";
  }
  
  
  print "$vcf max is $max_ID and $genesblock\n";
  #print "samtools faidx /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked $genesblock | bcftools consensus -H 1pIu $dir$a/$final_sample_name.vcf.gz > $output$a/$i/$final_sample_name.fasta\n";
  foreach $key (keys %phase_ID) {
    #print "$key : $phase_ID{$key} \n";
  }
  
  

} 

sub check_filter_2 {
  my ($vcf, $dir, $a, $genesblock, $filter) = @_;
  
  $genesblock =~ s/\R//g;
  #print "$vcf $genesblock\n";
  $total_loci = `bcftools view -H -r $genesblock $vcf | wc -l`;
  $miss = `bcftools view -H -r $genesblock -g miss $vcf| wc -l `;
  $total_loci =~ s/\R//g;
  $miss =~ s/\R//g;
  if(length($miss) == -1) {
    $miss = 0;
  }
  if(length($total_loci) == -1) {
    return "false";
  }
  #print "$vcf: total: $total_loci miss: $miss\n";
  $prop_miss = $miss/$total_loci * 100;
  #print "$vcf $prop_miss\n";
  if($prop_miss > 50) {
    return "false";
  }else{
    return "true";
  }

}