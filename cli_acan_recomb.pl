#!/usr/bin/perl
#counting allele freq of Clitarchus and Acanthoxyla in hybrid sample

use List::MoreUtils qw(uniq);


my ($vcf, $sample) = @ARGV;
$start = "false";


open(my $fh, '<:encoding(UTF-8)', $vcf)
  or die "Could not open file '$vcf' $!";
 
while (my $line = <$fh>) {
  $line =~ s/\R//g;
  @temp = split(/\s+/, $line);
  @cli = ();
  @acan = ();
  
  if($temp[0] eq "#CHROM") {
    for($i = 9; $i < $#temp + 1; $i ++) {
      #print "$i: $temp[$i]\n";
      if($temp[$i] eq "cli525" || $temp[$i] eq "cli753") {
        #print "cli at $i\n";
        push(@cli_sample, $i)
      }elsif($temp[$i] eq $sample) {
        #print "sample at $i\n";
        $sample_index = $i;
      }elsif($temp[$i] eq "axx" || $temp[$i] eq "axw2") {
        push(@acan_sample, $i)
      }
    }
  
    $start = "true";  
  }elsif($start eq "true") {
    
    foreach $a (@acan_sample) {
      push(@acan, get_allele($temp[$a]))
    }
    
    
    foreach $a (@cli_sample) {
      push(@cli, cli_allele($temp[$a]))
    }
   
    
    @unique_cli = uniq @cli;
    
    #print "target: at $sample_index $temp[$sample_index]\n";
    if($#unique_cli > -1) {
      plot($temp[$sample_index]);
    }
    
    
  
  }
  
  
  
}



sub get_allele() {
 
  my @locus = split(/:/, $_[0]);
  @result = ();
  if(index($locus[0], "/") != -1){
    @hap = split(/\//, $locus[0]);
    #print "$hap[0], $hap[1]\n";
  }else{
    @hap = split(/\|/, $locus[0]);
  }
  
  foreach $a (@hap) {
    if($a ne ".") {
      push(@result, $a);
    }
  
  }
  
  return @result;
}

sub cli_allele() {
 
  my @locus = split(/:/, $_[0]);
  @result = ();
  if(index($locus[0], "/") != -1){
    @hap = split(/\//, $locus[0]);
    #print "$hap[0], $hap[1]\n";
  }else{
    @hap = split(/\|/, $locus[0]);
  }
  
  foreach $a (@hap) {
    if($a ne ".") {
      my $pure_acan = "true";
      foreach $b (@acan) {
         
        if($a eq $b) {
          $pure_acan = "false";
          #print "untrue as $a and $b are the same\n";
          last;
        }
      }
      
      if($pure_acan eq "true") {
        push(@result, $a);
        #print "got in\n";
      }
    }
  
  }
  
  return @result;
}



sub plot() {

  my @locus = split(/:/, $_[0]);
  my $cli_freq = 0;
  my @recomb = ();
  my @non_recomb = ();
  my @non_missing = ();
  if(index($locus[0], "/") != -1){
    @hap = split(/\//, $locus[0]);
    #print "$hap[0], $hap[1]\n";
  }else{
      @hap = split(/\|/, $locus[0]);
  }

  foreach $a (@hap) {
    if($a ne ".") {
      #print "$a\n";
      push(@non_missing, $a);
    }
  
  }
  
  @unique_tar = uniq @non_missing;
  @hap_AD = split(/\,/, $locus[1]);

  $DP = $locus[2];
  
  my $cli_present = "false";
  #print $unique_tar;
  #print "\n";
  #print "$#unique_tar\n";
  if($#unique_tar > -1 && $DP > 0){
    
     
    foreach  $b (@unique_tar) {
      foreach $c (@unique_cli) {
       
        if($c eq $b) {
          $cli_present = "true";
          my $temp_freq = $hap_AD[$c]/$DP;
          $cli_freq = $cli_freq + $temp_freq;
          last;
       
        }
      }
    }
    
    if($sample eq "axx500") {
      
      if($cli_freq >= 0.35 && $cli_freq <= 0.65){
        print "recomb\tnorm\t$DP\n"  
      }elsif($cli_freq < 0.35) {
        print "non_recomb\tlow\t$DP\n";
      }elsif($cli_freq > 0.65) {
        print "non_recomb\thigh\t$DP\n";
      }
      
    }else{

      if ($cli_freq >= 0.1 && $cli_freq <= 0.5) {
        print "recomb\tnorm\t$DP\n"
      }elsif($cli_freq < 0.1) {
        print "non_recomb\tlow\t$DP\n";
      }elsif($cli_freq > 0.5) {
        print "non_recomb\thigh\t$DP\n";
      }
    }
  }
  


}


