#!/usr/bin/perl
#gettting heatmap? for database

@data = `cat $ARGV[0]`;

@pop = ("inermis_1.HEB", "inermis_2.HEB", "prasina_1.HEB", "prasina_2.HEB");
print "Genes\tA.inermis_1\tA.inermis_2\tA.prasina_1\tA.prasina_2\n";

foreach $a (@data) {
  $a =~ s/\R//g;
  print "$a\t";
  foreach $b (@pop) {
    $result = value($a, $b);
    print "$result\t";
  }
  print "\n";
} 


sub value () {
  #@genome = ("/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_acan/new_2/cor/");
  @genome = ("/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_acan/new_2/cor/", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_cli/new_2/cor/", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/arg100/HEB_2/cor/");
  
  $final_FPKM = 0;
  
  
  foreach $c (@genome) {
    #print "$c$_[1]\n";
    $temp_string = `grep -w "$_[0]" $c$_[1]`;
    $temp_FPKM = get_values($temp_string);
    
    $final_FPKM = $final_FPKM + $temp_FPKM;
  }
  
  return $final_FPKM;

}

sub get_values() {
  
  $_[0] =~ s/\R//g;
  @temp = split(/\s+/, $_[0]);
  
  
  if($temp[4] eq "hookeri") {
    $result = $temp[1];
  
  }elsif($temp[4] eq "geisovii") {
    $result = -$temp[2];
  }else{
    $result = 0;
  }
  
    return $result;


}