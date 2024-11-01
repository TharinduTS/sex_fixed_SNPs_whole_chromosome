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
Then copy the r script compare_Fasta to To_download directory
```
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
library(seqinr)
library(reshape2)
library(ggplot2)


#************************************
#install.packages("RFLPtools")
library(RFLPtools)

#create empty df to save final summary
# created vector with 3 characters
columns= c("Region","male_fixed_SNPs","female_fixed_SNPs") 

# and nrow with 0
myData = data.frame(matrix(nrow = 0, ncol = length(columns))) 

# assign column names
colnames(myData) = columns

chr_numbers<-1:20

for (region in chr_numbers) {
  current_region<-paste(region,"mb_vcfs_alignment.fasta",sep = '')
  print(current_region)

my_fasta<-read.fasta(current_region)
my_seq_df<-as.data.frame(my_fasta)

#Drop any rows with 'n' in any sample

#my_seq_df<-my_seq_df[rowSums(my_seq_df == "n")==0, , drop = FALSE]

#add location as a column
my_seq_df_with_loc <- tibble::rownames_to_column(my_seq_df, "location")

#export sequence from a selected sample to blast

# ********EDIT SELECTED SAMPLE ACCORDINGLY******

selected_sample<-'F_Ghana_WZ_BJE4687_combined__sorted.bam'

#**************************************************

full_sample_sequence<-my_seq_df_with_loc[[selected_sample]]

full_sequence<-paste(full_sample_sequence,collapse = '')

# remove markings for pure sequence

unmarked_full_sequence<-stringr::str_replace_all(full_sequence, '\\*', '')
unmarked_full_sequence<-stringr::str_replace_all(unmarked_full_sequence, '\\(', '')
unmarked_full_sequence<-stringr::str_replace_all(unmarked_full_sequence, '\\)', '')

writeLines(unmarked_full_sequence,"./full_sequence_of_selected_sample")
# 
# female_fixed <- my_seq_df[my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_Nigeria_EUA0333_combined__sorted.bam 
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Nigeria_EUA0334_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Nigeria_EUA0335_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
#                           & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]
# 
# male_fixed <- my_seq_df[my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_Nigeria_EUA0333_combined__sorted.bam 
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Nigeria_EUA0334_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Nigeria_EUA0335_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
#                         & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]


# Without Nigeria

female_fixed_no_niger <- my_seq_df[my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
                          #lab and JBL
                          #& my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$XT11_WW_trim_no_adapt_scafconcat_sorted.bam
                          #& my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$XT10_WZ_no_adapt._sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$JBL052_concatscafs_sorted.bam
                          # *******************
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
                          #lab and JBL
                          #& my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$XT7_WY_no_adapt__sorted.bam
                          #& my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$XT1_ZY_no_adapt._sorted.bam
                          # ******************
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]

#add a dummy row to avoid conflict if there are no female fixed SNPs

number_to_decuct_females<-0
if (nrow(female_fixed_no_niger)==0) {
  print(paste('No female fixed SNPs on',current_region))
  col_count<-ncol(female_fixed_no_niger_with_loc)-1
  empty_fix_row<-rep("X", each=col_count)
  number_to_decuct_females<-1
  
  female_fixed_no_niger[nrow(female_fixed_no_niger) + 1,]<-empty_fix_row
}

male_fixed_no_niger <- my_seq_df[my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam
                        #lab and JBL
                        #& my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$XT11_WW_trim_no_adapt_scafconcat_sorted.bam
                        #& my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$XT10_WZ_no_adapt._sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$JBL052_concatscafs_sorted.bam
                        # *******************
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
                        #lab and JBL
                        #& my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$XT7_WY_no_adapt__sorted.bam
                        #& my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$XT1_ZY_no_adapt._sorted.bam
                        # ******************
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]

#add a dummy row to avoid conflict if there are no male fixed SNPs

number_to_decuct_males<-0
if (nrow(male_fixed_no_niger)==0) {
  print(paste('No male fixed SNPs on',current_region))
  col_count<-ncol(male_fixed_no_niger_with_loc)-1
  empty_fix_row<-rep("X", each=col_count)
  number_to_decuct_males<-1
  
  male_fixed_no_niger[nrow(male_fixed_no_niger) + 1,]<-empty_fix_row
}


#mark locations with male and female fixed SNPs

library(tibble)
male_fixed_no_niger_with_loc <- tibble::rownames_to_column(male_fixed_no_niger, "location")

female_fixed_no_niger_with_loc <- tibble::rownames_to_column(female_fixed_no_niger, "location")


male_location_list<-male_fixed_no_niger_with_loc$location

female_location_list<-female_fixed_no_niger_with_loc$location

#Marking male fixed SNPs

for (i in 1:length(male_location_list)) {
  for (j in 1:length(my_seq_df_with_loc$location)) {
    if (my_seq_df_with_loc$location[j]==male_location_list[i]) {
      my_seq_df_with_loc[[selected_sample]][j]<-paste("*",my_seq_df_with_loc[[selected_sample]][j],"*",sep = "")
    }
  }
}


#Marking female fixed SNPs

for (i in 1:length(female_location_list)) {
  for (j in 1:length(my_seq_df_with_loc$location)) {
    if (my_seq_df_with_loc$location[j]==female_location_list[i]) {
      my_seq_df_with_loc[[selected_sample]][j]<-paste("(",my_seq_df_with_loc[[selected_sample]][j],")",sep = "")
    }
  }
}


# This region is for extracting sequences as 1000 bp regions

# Creating directory to save text files

dir.create('marked_sequence_files')
dir.create('unmarked_sequence_files')
#dir.create('to_blast')

locations_left<-as.numeric(male_location_list)
while (0<length(locations_left)) {
  print(paste("extracting region "))


SNP_location<-min(as.numeric(locations_left))

# Change these to change sequence length******
upstream<-100
downstream<-900

#********************

extract_begin<-SNP_location-upstream
extract_end<-SNP_location+downstream

extracted_sequence<-my_seq_df_with_loc[[selected_sample]][extract_begin:extract_end]

selected_sequence<-paste(extracted_sequence,collapse = '')

selected_sequence

print(paste("extracted SNPs from",extract_begin,"to",extract_end))

saving_file_name<-paste('Sex_sp_locations_marked_sequence_from_',extract_begin,'_to_',extract_end,'.txt',sep = '')
saving_file_name_unmarked<-paste('Sex_sp_locations_unmarked_sequence_from_',extract_begin,'_to_',extract_end,'.txt',sep = '')

print(paste("saving the sequence as",saving_file_name,"Males fixed SNPs are marked by * and Female fixed SNPs are marked by () in",selected_sample))

#removing file if already exists to avoid conflicts

file.remove(paste('./marked_sequence_files/',saving_file_name,sep = ''))
file.remove(paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''))
#file.remove(paste('./to_blast/',saving_file_name_unmarked,sep = ''))

#writing marked file only for selected file
first_line<-paste(">",selected_sample,sep = '')
write(first_line,file = paste('./marked_sequence_files/',saving_file_name,sep = ''),append = TRUE)
write(selected_sequence,file = paste('./marked_sequence_files/',saving_file_name,sep = ''),append = TRUE)

#writing unmarked file only for blast

# remove markings for pure sequence

# unmarked_selected_sequence<-stringr::str_replace_all(selected_sequence, '\\*', '')
# unmarked_selected_sequence<-stringr::str_replace_all(unmarked_selected_sequence, '\\(', '')
# unmarked_selected_sequence<-stringr::str_replace_all(unmarked_selected_sequence, '\\)', '')
# 
# write(first_line,file = paste('./to_blast/',saving_file_name_unmarked,sep = ''),append = TRUE)
# write(unmarked_selected_sequence,file = paste('./to_blast/',saving_file_name_unmarked,sep = ''),append = TRUE)

#writing unmarked multifasta for all samples

sample_list<-colnames(my_seq_df_with_loc)
sample_list<-sample_list[2:length(sample_list)]
for (sample in sample_list) {
  current_selected_sample<-sample
  first_line<-paste(">",current_selected_sample,sep = '')
  current_sequence<-my_seq_df_with_loc[[current_selected_sample]][extract_begin:extract_end]
  current_sequence_to_write<-paste(current_sequence,collapse = '')
  
  
  #for unmarked sequence file
  # remove markings for pure sequence
  
  unmarked_current_sequence_to_write<-stringr::str_replace_all(current_sequence_to_write, '\\*', '')
  unmarked_current_sequence_to_write<-stringr::str_replace_all(unmarked_current_sequence_to_write, '\\(', '')
  unmarked_current_sequence_to_write<-stringr::str_replace_all(unmarked_current_sequence_to_write, '\\)', '')
  
  #write unmarked files
  write(first_line,file = paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''),append = TRUE)
  write(unmarked_current_sequence_to_write,file = paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''),append = TRUE)
}

locations_left<-subset(locations_left,locations_left>extract_end)

}

#cleaning data

cleaned_male_fixed_no_niger_with_loc<-data.frame(male_fixed_no_niger_with_loc$location,male_fixed_no_niger_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam,male_fixed_no_niger_with_loc$F_SierraLeone_AMNH17272_combined__sorted.bam,male_fixed_no_niger_with_loc$F_SierraLeone_AMNH17274_combined__sorted.bam,male_fixed_no_niger_with_loc$F_IvoryCoast_xen228_combined__sorted.bam,male_fixed_no_niger_with_loc$JBL052_concatscafs_sorted.bam,male_fixed_no_niger_with_loc$M_Ghana_WY_BJE4362_combined__sorted.bam,male_fixed_no_niger_with_loc$M_Ghana_ZY_BJE4360_combined__sorted.bam,male_fixed_no_niger_with_loc$M_SierraLeone_AMNH17271_combined__sorted.bam,male_fixed_no_niger_with_loc$M_SierraLeone_AMNH17273_combined__sorted.bam)
colnames(cleaned_male_fixed_no_niger_with_loc)<-c("loc","F_Ghana_WZ_BJE4687_combined__sorted.bam","F_SierraLeone_AMNH17272_combined__sorted.bam","F_SierraLeone_AMNH17274_combined__sorted.bam","F_IvoryCoast_xen228_combined__sorted.bam","JBL052_concatscafs_sorted.bam","M_Ghana_WY_BJE4362_combined__sorted.bam","M_Ghana_ZY_BJE4360_combined__sorted.bam","M_SierraLeone_AMNH17271_combined__sorted.bam","M_SierraLeone_AMNH17273_combined__sorted.bam")
cleaned_male_fixed_no_niger_with_loc<-cleaned_male_fixed_no_niger_with_loc[order(cleaned_male_fixed_no_niger_with_loc$loc),]
write.table(cleaned_male_fixed_no_niger_with_loc, file='./cleaned_male_fixed_no_niger_with_loc.txt', quote=FALSE, sep='\t', col.names = NA)

cleaned_female_fixed_no_niger_with_loc<-data.frame(female_fixed_no_niger_with_loc$location,female_fixed_no_niger_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam,female_fixed_no_niger_with_loc$F_SierraLeone_AMNH17272_combined__sorted.bam,female_fixed_no_niger_with_loc$F_SierraLeone_AMNH17274_combined__sorted.bam,female_fixed_no_niger_with_loc$F_IvoryCoast_xen228_combined__sorted.bam._sorted.bam,female_fixed_no_niger_with_loc$JBL052_concatscafs_sorted.bam,female_fixed_no_niger_with_loc$M_Ghana_WY_BJE4362_combined__sorted.bam,female_fixed_no_niger_with_loc$M_Ghana_ZY_BJE4360_combined__sorted.bam,female_fixed_no_niger_with_loc$M_SierraLeone_AMNH17271_combined__sorted.bam,female_fixed_no_niger_with_loc$M_SierraLeone_AMNH17273_combined__sorted.bam)
colnames(cleaned_female_fixed_no_niger_with_loc)<-c("loc","F_Ghana_WZ_BJE4687_combined__sorted.bam","F_SierraLeone_AMNH17272_combined__sorted.bam","F_SierraLeone_AMNH17274_combined__sorted.bam","F_IvoryCoast_xen228_combined__sorted.bam","JBL052_concatscafs_sorted.bam","M_Ghana_WY_BJE4362_combined__sorted.bam","M_Ghana_ZY_BJE4360_combined__sorted.bam","M_SierraLeone_AMNH17271_combined__sorted.bam","M_SierraLeone_AMNH17273_combined__sorted.bam")
cleaned_female_fixed_no_niger_with_loc<-cleaned_female_fixed_no_niger_with_loc[order(cleaned_female_fixed_no_niger_with_loc$loc),]
write.table(cleaned_female_fixed_no_niger_with_loc, file='./cleaned_female_fixed_no_niger_with_loc.txt', quote=FALSE, sep='\t', col.names = NA)

# without lab and JBL
# cleaned_female_fixed_no_niger_with_loc<-data.frame(female_fixed_no_niger_with_loc$location,female_fixed_no_niger_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam,female_fixed_no_niger_with_loc$F_SierraLeone_AMNH17272_combined__sorted.bam,female_fixed_no_niger_with_loc$F_SierraLeone_AMNH17274_combined__sorted.bam,female_fixed_no_niger_with_loc$F_IvoryCoast_xen228_combined__sorted.bam,female_fixed_no_niger_with_loc$M_Ghana_WY_BJE4362_combined__sorted.bam,female_fixed_no_niger_with_loc$M_Ghana_ZY_BJE4360_combined__sorted.bam,female_fixed_no_niger_with_loc$M_SierraLeone_AMNH17271_combined__sorted.bam,female_fixed_no_niger_with_loc$M_SierraLeone_AMNH17273_combined__sorted.bam)
# colnames(cleaned_female_fixed_no_niger_with_loc)<-c("loc","F_Ghana_WZ_BJE4687_combined__sorted.bam","F_SierraLeone_AMNH17272_combined__sorted.bam","F_SierraLeone_AMNH17274_combined__sorted.bam","F_IvoryCoast_xen228_combined__sorted.bam","M_Ghana_WY_BJE4362_combined__sorted.bam","M_Ghana_ZY_BJE4360_combined__sorted.bam","M_SierraLeone_AMNH17271_combined__sorted.bam","M_SierraLeone_AMNH17273_combined__sorted.bam")
# cleaned_female_fixed_no_niger_with_loc<-cleaned_female_fixed_no_niger_with_loc[order(cleaned_female_fixed_no_niger_with_loc$loc),]
# write.table(cleaned_female_fixed_no_niger_with_loc, file='./cleaned_female_fixed_no_niger_with_loc.txt', quote=FALSE, sep='\t', col.names = NA)

# make seperate directory for current region
region_name<-paste(current_region,'_summaries',sep = '')
dir.create(region_name)


#move summary files there
file.copy(from = "cleaned_female_fixed_no_niger_with_loc.txt",
          to   = region_name)
file.remove("cleaned_female_fixed_no_niger_with_loc.txt")

file.copy(from = "cleaned_male_fixed_no_niger_with_loc.txt",
          to   = region_name)
file.remove("cleaned_male_fixed_no_niger_with_loc.txt")

file.copy(from = "full_sequence_of_selected_sample",
          to   = region_name)
file.remove("full_sequence_of_selected_sample")

#move folders
file.rename('./marked_sequence_files',paste(region_name,"/marked_sequence_files",sep = ''))
file.rename('./unmarked_sequence_files',paste(region_name,"/unmarked_sequence_files",sep = ''))

#save details in the final summary file

region_save_name<-paste(region,sep = '')
male_fixed_number<-nrow(cleaned_male_fixed_no_niger_with_loc)-number_to_decuct_males
female_fixed_number<-nrow(cleaned_female_fixed_no_niger_with_loc)-number_to_decuct_females


myData[nrow(myData) + 1,]<-c(region_save_name,male_fixed_number,female_fixed_number)


}

#print final summary and save it
print(myData)
write.table(myData, file='./all_region_summary.txt', quote=FALSE, sep='\t', col.names = NA)

#plot results

# Reshape data from wide to long format
long_df <- melt(myData, id.vars = c("Region(mb)"), variable.name = "Sex", 
                value.name = "Fixed_SNPs")

long_df<-long_df[order(long_df$Region),]

fixed_snp_plot<-ggplot(long_df, aes(fill=Sex, y=Fixed_SNPs, x=Region)) +
  geom_bar(position='dodge', stat='identity')+
  theme_bw()+
  theme_classic()

ggsave("fixed_snp_plot.pdf",fixed_snp_plot)

```

Then make sure all the needed R packages are installed first by running R script and installing needed packages
```
install.packages("rstudioapi")
install.packages("seqinr")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("ggplot2")

install.packages("RFLPtools")

```
Then make new files for different chromosome regions
```
for i in {1..20};do sed "s/chr_numbers<-1/chr_numbers<-${i}/g" compare_Fasta>compare_Fasta${i}.sh;done
```
move those fioles into corresponding directories
```
for i in {1..20}; do mkdir region${i}; mv ${i}mb_vcfs_alignment.fasta region${i}; mv compare_Fasta${i}.sh region${i ;done
```
Change rscript directory accordingly
```
for i in {1..20}; do sed -i "s/testing_whole_chr7\/to_download/testing_whole_chr7\/to_download\/region${i}/g" region${i}/compare_Fasta${i}.sh;done
```
And run them
```
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=64gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben
#SBATCH --array=1-20

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load StdEnv/2020 r/4.3.1

Rscript region${SLURM_ARRAY_TASK_ID}/compare_Fasta${SLURM_ARRAY_TASK_ID}.sh
```

