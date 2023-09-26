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
  my $output = join "", $wrkdir, "/slurm/axw2_", $name, ".sl";
  my $dir = join "", $wrkdir,$i, "\/";
  my $scaff = join "", $temp[0], "\:", $temp[1],"-", $temp[2];
  $mini = "mini_";

 
  print "preparing gene$i\n";
 
  #`samtools faidx /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked $scaff \> $dir$name.fasta`;
  #`samtools faidx $dir$name.fasta`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/genome_axw2/all_AXW2_2.vcf.gz -r $scaff > $dir/AXW2_2.vcf`;
  



  unless(open FILE, '>'.$output) {
	  # Die with error message 
	  # if we can't open it.
  	die "nUnable to create $output/n";
  }

  print FILE "#!/bin/bash -e\n";
  print FILE "#SBATCH -J AXW2_$name\n";
  print FILE "#SBATCH -A ga02470\n";
  print FILE "#SBATCH --ntasks=1\n";
  print FILE "#SBATCH --cpus-per-task=1\n";
  print FILE "#SBATCH --mem=1G\n";
  print FILE "#SBATCH --time=10:00\n";
  print FILE "#SBATCH --hint=nomultithread\n";
  print FILE "#SBATCH --partition=large\n";
  print FILE "#SBATCH -o $output.txt\n";
  print FILE "#SBATCH -e $output.err\n";
  print FILE "#SBATCH --mail-type=FAIL\n";
  print FILE "#SBATCH --mail-user=scho249\@aucklanduni.ac.nz\n";
  print FILE "#SBATCH --chdir=$dir\n\n";
  
print FILE "ml purge\n";
print FILE "ml load Miniconda3/22.11.1-1\n";
print FILE "source activate whatshap-env\n\n";
  
print FILE "whatshap phase -o phased_AXW2_2.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked AXW2_2.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/genome_axw2/add_AXW2.bam --indels\n";

print FILE "whatshap stats --block-list=phased_AXW2_2.tsv --gtf=phased_AXW2_2.gtf phased_AXW2_2.vcf\n";


}

