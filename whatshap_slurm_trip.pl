#!/usr/bin/perl
#running sample to whatshap in parallel
#commandL perl "/scale_wlg_nobackup/filesets/nobackup/ga02470/prevent_from_delete/script/whatshap_slurm_trip.pl" "/scale_wlg_nobackup/filesets/nobackup/ga02470/freebayes/lowCov/freebayes_3/whatshap/busco/all/list.tsv" "/scale_wlg_nobackup/filesets/nobackup/ga02470/freebayes/lowCov/freebayes_3/whatshap/busco/all/slurm/" /scale_wlg_nobackup/filesets/nobackup/ga02470/freebayes/lowCov/freebayes_3/whatshap/busco/all/both_indel/


$wrkdir = $ARGV[0];

@genes = `cat "/scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/C.hooker_busco/run_insecta_odb10/single_region.bed"`;

for($i = 0; $i < $#genes + 1; $i ++) {
  `mkdir $wrkdir/$i/`;
  
  $genes[$i] =~ s/\R//g;
  my @temp = split(/\s+/, $genes[$i]);
  my $name = join "","gene", $i; 
  my $output = join "", $wrkdir, "/slurm/", $name, ".sl";
  my $dir = join "", $wrkdir,$i, "\/";
  my $scaff = join "", $temp[0], "\:", $temp[1],"-", $temp[2];
  $mini = "mini_";

 
  print "preparing gene$i\n";
  #`samtools faidx /scale_wlg_nobackup/filesets/nobackup/ga02470/freebayes/gapclosed.fasta.masked $scaff \> $dir$name.fasta`;
  #`samtools faidx $dir$name.fasta`;
  

=pod 
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_axx.vcf.gz -r $scaff > $dir/axx.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_axp.vcf.gz -r $scaff > $dir/axp.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_axg.vcf.gz -r $scaff > $dir/axg.vcf`;
=cut

  unless(open FILE, '>'.$output) {
	  # Die with error message 
	  # if we can't open it.
  	die "nUnable to create $result/n";
  }

  print FILE "#!/bin/bash -e\n";
  print FILE "#SBATCH -J $name\n";
  print FILE "#SBATCH -A ga02470\n";
  print FILE "#SBATCH --ntasks=1\n";
  print FILE "#SBATCH --cpus-per-task=30\n";
  print FILE "#SBATCH --mem=15G\n";
  print FILE "#SBATCH --time=72:00:00\n";
  print FILE "#SBATCH --hint=nomultithread\n";
  print FILE "#SBATCH -o $output.txt\n";
  print FILE "#SBATCH -e $output.err\n";
  print FILE "#SBATCH --mail-type=FAIL\n";
  print FILE "#SBATCH --mail-user=scho249\@aucklanduni.ac.nz\n";
  print FILE "#SBATCH --chdir=$dir\n\n";
  
  print FILE "ml purge\n";
  print FILE "ml load Miniconda3/22.11.1-1\n";
  print FILE "source activate whatshap-env\n\n";

  print FILE "whatshap polyphase -t 30 -p 3 -o phased_axg.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axg.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/ON/add_axg.bam --indels\nwhatshap stats --block-list=phased_axg.tsv --gtf=phased_axg.gtf phased_axg.vcf\n";

  print FILE "whatshap polyphase -t 30 -p 3 -o phased_axp.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axp.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/ON/add_axp.bam --indels\nwhatshap stats --block-list=phased_axp.tsv --gtf=phased_axp.gtf phased_axp.vcf\n";

  print FILE "whatshap polyphase -t 30 -p 3 -o phased_axx.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axx.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/ON/add_axx.bam --indels\nwhatshap stats --block-list=phased_axx.tsv --gtf=phased_axx.gtf phased_axx.vcf\n";
  

  
=pod  
  print FILE "whatshap polyphase -t 30 -p 3 -o phased_axg.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axg.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_axg.bam /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/ON/add_axg.bam --indels\nwhatshap stats --block-list=phased_axg.tsv --gtf=phased_axg.gtf phased_axg.vcf\n";

  print FILE "whatshap polyphase -t 30 -p 3 -o phased_axp.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axp.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_axp.bam /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/ON/add_axp.bam --indels\nwhatshap stats --block-list=phased_axp.tsv --gtf=phased_axp.gtf phased_axp.vcf\n";

  print FILE "whatshap polyphase -t 30 -p 3 -o phased_axx.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axx.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_axx.bam /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/ON/add_axx.bam --indels\nwhatshap stats --block-list=phased_axx.tsv --gtf=phased_axx.gtf phased_axx.vcf\n";
  
=cut  



}

