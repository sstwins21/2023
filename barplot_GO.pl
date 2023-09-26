#!/usr/bin/perl
#transform the output csv files from ShinyGo to format used in R for horizonatl barplot
use List::MoreUtils qw(uniq);

@data = `awk -F "," '{if(\$1 <= 0.05 && index(\$6, "amigo")) print \$6}' acan_inermis_1.csv  acan_inermis_2.csv  acan_prasina_1.csv  acan_prasina_2.csv  cli_inermis_1.csv  cli_inermis_2.csv  cli_prasina_1.csv  cli_prasina_2.csv | sort | uniq -c | awk '{if(\$1 >1) print \$2}' | awk -F "/" '{print \$NF}'`;
#@data = `awk -F "," '{if(\$1 <= 0.05 && index(\$6, "amigo")) print \$5}' acan_inermis_1.csv  acan_inermis_2.csv  acan_prasina_1.csv  acan_prasina_2.csv  cli_inermis_1.csv  cli_inermis_2.csv  cli_prasina_1.csv  cli_prasina_2.csv | sort | uniq -d`;
@pathways = ();


@pop =  ("/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_inermis_1.csv", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_inermis_2.csv", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_prasina_1.csv", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_prasina_2.csv");

#print "$pop[0]\n";


print "GO,A. inermis_1,A. inermis_2, A. prasina_1, A. prasina_2\n";
foreach $a (@data) {
  chomp $a;
  @temp  = split(/\s+/, $a);
  print "$a";
    
  foreach $b (@pop) {
    $fold = `grep "$a" $b | awk -F "," '{print \$4}'`;
    chomp $fold;
  
    if($fold ne "") {
      print ",$fold";
    }else{
      print ",0";
    }
  }

  print "\n";
}
