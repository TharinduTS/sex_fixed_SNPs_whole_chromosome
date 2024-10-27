# sex_fixed_SNPs_whole_chromosome


 Working in beluga /scratch/premacht/testing_whole_chr7

Make bed files for each mb in Chr 7
```
for i in {1..20} ; do echo "Chr7 ${i}000000 ${i}999999" >>${i}mb_region.bed;done
```

Then extracted regions

First I saved the sample extracting script as “extracting_regions.sh”


```
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load gatk/4.4.0.0 StdEnv/2023

for i in ../new_project_Apr_2023/new_data_Jun23/*.bam; do gatk HaplotypeCaller -R ../new_project_Apr_2023/new_data_Jun23/reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -I ${i} -L YOURNOmb_region.bed -ERC BP_RESOLUTION -O ${i##../new_project_Apr_2023/new_data_Jun23/}YOURNOmb.vcf;done

mkdir YOURNOmb_vcfs
mv *YOURNOmb.vcf* 1mb_vcfs
```

Then I used following command to make copies for each mb

Use double quotes here to let shell change the value for $i
```
for i in {1..20};do sed "s/YOURNO/$i/g" extracting_regions.sh >extract${i}mb.sh;done
```
Then submit all jobs
```
for i in {1..20};do sbatch extract${i}mb.sh ;done
```
First copy vcf2phylip to the main directory from 

https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py

and then run following to copy the code to each of the sub directories 
```
```
for i in 1; do cd ${i}mb_vcfs;cp ../vcf2phylip.py .;cd .. ;done
```
and then copy this script to the main directory as convert_vcf_and_align_fasta.sh
```
#!/bin/sh
#SBATCH --job-name=align_fasta
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=32gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

for i in *vcf;do python ../vcf2phylip.py --input ${i} --fasta;done
cat *.fasta >${PWD##*/}_combined_fasta.fa

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 mafft-mpi/7.471
mafft --auto ${PWD##*/}_combined_fasta.fa >${PWD##*/}_alignment.fasta
```
Then this code copies this to all the sub directories and submits jobs

```
for i in {1..20}; do cd ${i}mb_vcfs;cp ../convert_vcf_and_align_fasta.sh .;sbatch convert_vcf_and_align_fasta.sh;cd .. ;done
```


Then copied all the aligned files to a separate directory to download
```
mkdir to_download
for i in {1..20};do cp ${i}mb_vcfs/*alignment.fasta ./to_download/;done
```
Then download 


if (nrow(male_fixed_no_niger)==0) {
  print(paste('No male fixed SNPs on',current_region))
  col_count<-ncol(male_fixed_no_niger_with_loc)-1
  empty_fix_row<-rep("X", each=col_count)
  
  male_fixed_no_niger[nrow(male_fixed_no_niger) + 1,]<-empty_fix_row
}
```
