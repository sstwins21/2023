#!/usr/bin/perl
#getting gene IDs


@data = `cat $ARGV[0]`;

for $a (@data) {
  $a =~ s/\R//g;
  @temp = split(/\s+/, $a);
  if($ARGV[1] eq "arg100") {
    $result = `grep -w $temp[0] "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/arg100_david_2.txt"`;
  }elsif($ARGV[1] eq "acan") {
    $result = `grep -w $temp[0] "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/acan_david_2.txt"`;
  }elsif($ARGV[1] eq "cli") {
    $result = `grep -w $temp[0] "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/berkley/david_2.txt"`;
}
  @temp2 = split(/\s+/, $result);
  
  if($temp2[0] ne "") {
    $temp[1] =~ s/\R//g;
    if($ARGV[1] eq "gene") {
      print "$temp[0]\t$temp2[1]\n";
    }else{
      print "$temp2[1]\n";
    }
  }

}