#!/usr/bin/perl
#producing assembly qual slurm files 

#@sample = `cat "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/samples.txt"`;
@sample = ("axg", "axp", "axx500");


for($i = 1; $i < 19; $i ++) {
  slurm($i);
};


 


sub slurm () {

  my $num = $_[0];
  #my $output = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/window_slide/", $name,"/ksmooth/", $name, "_", $num, ".R";
  my $output = "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/window_slide/all/all_$num.R";
  my $dir = join "", "/scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/window_slide/", $name,"/ksmooth/";
  my $data = join "", $name, "_ad.txt";
  my $scaff = join "", "scaffold_", $num;
  my $plot_name = join "", "scaffold_", $num, ".pdf";

  unless(open FILE, '>'.$output) {
	  # Die with error message 
	  # if we can't open it.
  	die "nUnable to create $name/n";
  } 
  
  print FILE "library(ggplot2)\n";
  print FILE "library(KernSmooth)\n";
  print FILE "library(ggtext)\n\n";
  print FILE "axg_data <- read.table(\"axg_ad.txt\")\n";
  print FILE "axp_data <- read.table(\"axp_ad.txt\")\n";
  print FILE "axx500_data <- read.table(\"axx500_ad.txt\")\n";

  print FILE "axg_scaf <- axg_data[ which(axg_data\$V1==\'$scaff\'), ]\n";
  print FILE "axg_band <- dpik(axg_scaf\$V2, scalest = \"minim\", level = 2L, kernel = \"normal\",  canonical = FALSE, gridsize = 401L)\n";
  print FILE "axg_ksmooth_scaf <- ksmooth(axg_scaf\$V2, axg_scaf\$V3, kernel = \"normal\", bandwidth = axg_band)\n";

  print FILE "axg_newDT <- data.frame(axg_ksmooth_scaf)\n";
  print FILE "colnames(axg_newDT) <- c(\"position\", \"C.hookeri\")\n";
  
  print FILE "axp_scaf <- axp_data[ which(axp_data\$V1==\'$scaff\'), ]\n";
  print FILE "axp_band <- dpik(axp_scaf\$V2, scalest = \"minim\", level = 2L, kernel = \"normal\",  canonical = FALSE, gridsize = 401L)\n";
  print FILE "axp_ksmooth_scaf <- ksmooth(axp_scaf\$V2, axp_scaf\$V3, kernel = \"normal\", bandwidth = axp_band)\n";

  print FILE "axp_newDT <- data.frame(axp_ksmooth_scaf)\n";
  print FILE "colnames(axp_newDT) <- c(\"position\", \"C.hookeri\")\n";
  
  print FILE "axx500_scaf <- axx500_data[ which(axx500_data\$V1==\'$scaff\'), ]\n";
  print FILE "axx500_band <- dpik(axx500_scaf\$V2, scalest = \"minim\", level = 2L, kernel = \"normal\",  canonical = FALSE, gridsize = 401L)\n";
  print FILE "axx500_ksmooth_scaf <- ksmooth(axx500_scaf\$V2, axx500_scaf\$V3, kernel = \"normal\", bandwidth = axx500_band)\n";

  print FILE "axx500_newDT <- data.frame(axx500_ksmooth_scaf)\n";
  print FILE "colnames(axx500_newDT) <- c(\"position\", \"C.hookeri\")\n\n";
  
  print FILE "pdf(\"axg_scaffold_$num.pdf\", width=12)\n";
  print FILE "ggplot(data = axg_newDT, mapping = aes(x = position))  +  geom_line(aes(y = C.hookeri)) +  ylim(0, 1) + theme_classic() + geom_hline(yintercept=0.33, linetype=\"dashed\", color = \"red\") + labs(title= \'*Acanthoxyla prasina* (AXG)\', x = \"Position along scaffold (mbp)\", y = \'Frequency of *C. hookeri* diagnostic alleles\') + theme( plot.margin = margin(, 2, , , \"cm\"), text = element_text(size = 20), plot.title = element_markdown(),  axis.title.y = element_markdown()) + scale_x_continuous(labels = scales::label_number(suffix = \" mbp\", scale = 1e-6))\n";
  print FILE "dev.off()\n\n";
  
  print FILE "pdf(\"axp_scaffold_$num.pdf\", width=12)\n";
  print FILE "ggplot(data = axp_newDT, mapping = aes(x = position))  +  geom_line(aes(y = C.hookeri)) +  ylim(0, 1) + theme_classic() + geom_hline(yintercept=0.33, linetype=\"dashed\", color = \"red\") + labs(title= \'*Acanthoxyla prasina* (AXP)\', x = \"Position along scaffold (mbp)\", y = \'Frequency of *C. hookeri* diagnostic alleles\') + theme( plot.margin = margin(, 2, , , \"cm\"), text = element_text(size = 20), plot.title = element_markdown(),  axis.title.y = element_markdown()) + scale_x_continuous(labels = scales::label_number(suffix = \" mbp\", scale = 1e-6))\n";
  print FILE "dev.off()\n\n";
  
  print FILE "pdf(\"axx500_scaffold_$num.pdf\", width=12)\n";
  print FILE "ggplot(data = axx500_newDT, mapping = aes(x = position))  +  geom_line(aes(y = C.hookeri)) +  ylim(0, 1) + theme_classic() + geom_hline(yintercept=0.5, linetype=\"dashed\", color = \"red\") + labs(title= \'*Acanthoxyla inermis*\', x = \"Position along scaffold (mbp)\", y = \'Frequency of *C. hookeri* diagnostic alleles\') + theme( plot.margin = margin(, 2, , , \"cm\"), text = element_text(size = 20), plot.title = element_markdown(),  axis.title.y = element_markdown()) + scale_x_continuous(labels = scales::label_number(suffix = \" mbp\", scale = 1e-6))\n";
  print FILE "dev.off()\n";

  
}