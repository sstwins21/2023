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
  my $output = join "", $wrkdir, "/slurm/dip_", $name, ".sl";
  my $dir = join "", $wrkdir,$i, "\/";
  my $scaff = join "", $temp[0], "\:", $temp[1],"-", $temp[2];
  $mini = "mini_";

 
  print "preparing gene$i\n";
 
  `samtools faidx /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked $scaff \> $dir$name.fasta`;
  `samtools faidx $dir$name.fasta`;
  

=pod
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_axx500.vcf.gz -r $scaff > $dir/axx500.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_arg100.vcf.gz -r $scaff > $dir/arg100.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_axw2.vcf.gz -r $scaff > $dir/axw2.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_cli525.vcf.gz -r $scaff > $dir/cli525.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_cli753.vcf.gz -r $scaff > $dir/cli753.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_cli948.vcf.gz -r $scaff > $dir/cli948.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_pki14.vcf.gz -r $scaff > $dir/pki14.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_pss.vcf.gz -r $scaff > $dir/pss.vcf`;
  `bcftools view /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/vcfs/fixed_tp.vcf.gz -r $scaff > $dir/tp.vcf`;





  unless(open FILE, '>'.$output) {
	  # Die with error message 
	  # if we can't open it.
  	die "nUnable to create $result/n";
  }

  print FILE "#!/bin/bash -e\n";
  print FILE "#SBATCH -J raw_dip$name\n";
  print FILE "#SBATCH -A ga02470\n";
  print FILE "#SBATCH --ntasks=1\n";
  print FILE "#SBATCH --cpus-per-task=1\n";
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
  
print FILE "whatshap phase -o phased_axx500.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axx500.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_axx500.bam --indels\n";
print FILE "whatshap phase -o phased_arg100.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked arg100.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_arg100.bam --indels\n";
print FILE "whatshap phase -o phased_axw2.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked axw2.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_axw2.bam --indels\n";
print FILE "whatshap phase -o phased_cli525.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked cli525.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_cli525.bam --indels\n";
print FILE "whatshap phase -o phased_cli753.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked cli753.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_cli753.bam --indels\n";
print FILE "whatshap phase -o phased_cli948.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked cli948.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/uniq_cli948.bam --indels\n";
print FILE "whatshap phase -o phased_pki14.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked pki14.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_pki14.bam --indels\n";
print FILE "whatshap phase -o phased_pss.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked pss.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/add_pss.bam --indels\n";
print FILE "whatshap phase -o phased_tp.vcf --reference /scale_wlg_nobackup/filesets/nobackup/ga02470/annotation/gapclosed.fasta.masked tp.vcf /scale_wlg_nobackup/filesets/nobackup/ga02470/article/chapter3/map/uniq_tp.bam --indels\n";
print FILE "whatshap stats --block-list=phased_axx500.tsv --gtf=phased_axx500.gtf phased_axx500.vcf\n";
print FILE "whatshap stats --block-list=phased_arg100.tsv --gtf=phased_arg100.gtf phased_arg100.vcf\n";
print FILE "whatshap stats --block-list=phased_axw2.tsv --gtf=phased_axw2.gtf phased_axw2.vcf\n";
print FILE "whatshap stats --block-list=phased_cli525.tsv --gtf=phased_cli525.gtf phased_cli525.vcf\n";
print FILE "whatshap stats --block-list=phased_cli753.tsv --gtf=phased_cli753.gtf phased_cli753.vcf\n";
print FILE "whatshap stats --block-list=phased_cli948.tsv --gtf=phased_cli948.gtf phased_cli948.vcf\n";
print FILE "whatshap stats --block-list=phased_pki14.tsv --gtf=phased_pki14.gtf phased_pki14.vcf\n";
print FILE "whatshap stats --block-list=phased_pss.tsv --gtf=phased_pss.gtf phased_pss.vcf\n";
print FILE "whatshap stats --block-list=phased_tp.tsv --gtf=phased_tp.gtf phased_tp.vcf\n";
=cut
}

