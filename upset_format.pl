#!/usr/bin/perl
#transform the output csv files from ShinyGo to format used in R for upset plot
use List::MoreUtils qw(uniq);

@data = `awk -F "," '{if(\$1 <= 0.05 && index(\$6, "amigo")) print \$5}' acan_inermis_1.csv  acan_inermis_2.csv  acan_prasina_1.csv  acan_prasina_2.csv  cli_inermis_1.csv  cli_inermis_2.csv  cli_prasina_1.csv  cli_prasina_2.csv | sort | uniq -c | sort`;
@pathways = ();


@pop =  ("/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_inermis_1.csv", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_inermis_2.csv", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_prasina_1.csv", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter4/GO/" . $ARGV[0] . "_prasina_2.csv");

#print "$pop[0]\n";


print "Pathways,A. inermis_1,A. inermis_2, A. prasina_1, A. prasina_2\n";
foreach $a (@data) {
  chomp $a;
  @temp  = split(/\s+/, $a);
  
  if($temp[1] > 1) {
    $pathways = "";
    for($i = 2; $i < $#temp + 1; $i ++) {
      $pathways = "$pathways$temp[$i] ";
    }
    $pathways =~ s/\s+$//;
    push(@pathways, @pathways);
    #print "$pathways,";
    for($i = 0; $i < $#pop + 1; $i ++) {
      if($i > 0) {
        print ",";
      }
      $geneCount = `awk -F "," '{if(\$1 <= 0.05 &&\$5 == "$pathways" && index(\$6, "amigo")) print \$4}' $pop[$i]`;
      $GO = `awk -F "," '{if(\$5 == "$pathways" && index(\$6, "amigo")) print \$6}' $pop[$i]`;
      if($geneCount ne "") {
        chomp $geneCount;
        #print "$geneCount";
        print "1"
      }else{
        print "0";
      }
    }
    print"\n";
  }
  
}
