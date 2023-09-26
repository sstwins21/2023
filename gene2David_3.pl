#!/usr/bin/perl
#getting HEB for single genes



($data, $genome, $pop) = @ARGV;

@data = `cat $data`;
%HEB = ();
$gei = "geisovii";
$hoo = "hookeri";
$bal = "balanced";

if($pop eq "inermis_1") {
  $pop_index = 3;
}elsif($pop eq "inermis_2") {
  $pop_index = 4;
}elsif($pop eq "prasina_1") {
  $pop_index = 5;
}elsif($pop eq "prasina_2") {
  $pop_index = 6;
}

for $a (@data) {
  $a =~ s/\R//g;
  @temp = split(/\s+/, $a);
  if($ARGV[1] eq "arg100") {
    $result = `grep -w $temp[0] "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/arg100_david.txt"`;
    $pop_count = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/arg100/HEB/arg100_FPKM.txt";
    
  }elsif($ARGV[1] eq "acan") {
    $result = `grep -w $temp[0] "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/acan_david.txt"`;
    $pop_count = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_acan/new/acan_FPKM.txt"
    
  }elsif($ARGV[1] eq "cli") {
    $result = `grep -w $temp[0] "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/david.txt"`;
    $pop_count = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/wasp_cli/new/cli_FPKM.txt";
  }  
  
  $result =~ s/\R//g;
  @temp2 = split(/\s+/, $result);
  if($temp2[1] ne "") {
    if (!exists $HEB{$temp2[1]}) {
      $HEB{$temp2[1]}{"geisovii"} = 0;
      $HEB{$temp2[1]}{"hookeri"} = 0;
      $HEB{$temp2[1]}{"balanced"} = 0;
    }
  
    $count = `grep -w $temp2[0] $pop_count`;
    $count =~ s/\R//g;
    @temp3 = split(/\s+/, $count);
    
    $HEB{$temp2[1]}{$temp[1]} = $HEB{$temp2[1]}{$temp[1]} + ($temp3[$pop_index] * $temp[5]);
    #print "$result $temp[1]\n";
    #print "$HEB{$temp2[1]}{$temp[1]} $HEB{$temp2[1]}{$gei} $temp[1]\n";
    #print "$count num: $pop_index $pop index: $temp3[$pop_index]\n";
  }
}

foreach my $key (keys %HEB) {

  if($HEB{$key}{$hoo} > 0) {
    
    if($HEB{$key}{$gei} > 0) {
      
      if($HEB{$key}{$hoo}/(($HEB{$key}{$hoo} + $HEB{$key}{$gei})) > 0.66) {
        $final_HEB = "hookeri";
        #$final_ratio = $HEB{$key}{$hoo}/($HEB{$key}{$hoo} + $HEB{$key}{$gei})
      }elsif($HEB{$key}{$gei}/(($HEB{$key}{$hoo} + $HEB{$key}{$gei})) > 0.66) {
        $final_HEB = "geisovii";
        #$final_ratio = $HEB{$key}{$gei}/($HEB{$key}{$hoo} + $HEB{$key}{$gei})
      }else{
        $final_HEB = "balanced";
      }
      
      
    }else{
      $final_HEB = "hookeri";
    }
  
  }elsif($HEB{$key}{$gei} > 0) {
    $final_HEB = "geisovii";
  }else{
    $final_HEB = "balanced";
  }

  print "$key\t$HEB{$key}{$hoo}\t$HEB{$key}{$gei}\t$HEB{$key}{$bal}\t$final_HEB\n"; 
}

