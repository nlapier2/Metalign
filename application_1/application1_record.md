# Application 1 record

### Contents:

- [Analysis logic](#analysis)
- [Data exploration](#explore)
- [Process all data](#runall)
- [Summarizing results](#summary)



---

### Analysis logic: <a name="analysis"></a>

1. Data exploration
   - [x] Compare the 2 rep of "Pig_14_Day_14_R1" data
   - [x] Compare the merged R1R2 of the 2 rep of "Pig_14_Day_14" data
   - [ ] Pick the top extra species in merged_R1R2_combine data, align the 3 reads file to the genome 
2. Run all data
   - [ ] 
3. Summarizing results

---

### Data exploration: <a name="explore"></a>

1. Compare 2 read1 files of "Pig_14_Day_14_R1"

   - code:

   ```bash
   #script: a1.3.2_run_combined_keepdup.pbs
   metalign_py="/gpfs/scratch/sml6467/tools/Metalign/metalign.py"
   metalign_data="/storage/home/sml6467/scratch/tools/Metalign/data"
   python3 $metalign_py <fastq_file>  $metalign_data  --output <output_name> --keep_temp_files --temp_dir <output_temp_folder>
   ```

   - results:
     - [original vs reseq](https://drive.google.com/file/d/1DUzZwB4OxxFc2S-g0H4vXP6EXBrRfMVb/view)
     - [combined vs original+reseq](https://drive.google.com/file/d/1OPzYJWfABIjTKW0zB0ckLJ-ulRhkbc7b/view)

   - conclusion: good rep consistency, but combined data has more species



2. Compare merge-R1R2 data of "Pig_14_Day_14_R1"

   - code

   ```bash
   # merge all the R1 R2 data (code for combined folder, rest are same)
   for file in `ls *R1*`; do
       d_pre=`echo ${file%_R1*}`
       d_suf=`echo ${file#*R1}`
       echo ${d_pre}_R1${d_suf} ${d_pre}_R2${d_suf}
       cat ${d_pre}_R1${d_suf} ${d_pre}_R2${d_suf} > merged_R12_${d_pre}_combined${d_suf}
       unset file d_pre d_suf
   done
   
   # run metalign 
   #script: a1.5_run_merged_R12_files.pbs
   metalign_py="/gpfs/scratch/sml6467/tools/Metalign/metalign.py"
   metalign_data="/storage/home/sml6467/scratch/tools/Metalign/data"
   cd /gpfs/group/dmk333/default/shaopeng/projects/202002_metalign_application/results/20200215_application1_variation/merged_R1R2_output
   
   for file in `ls *.fastq`; do
     echo "Processing $file"
     name=`echo ${file##*_}`
     name=`echo ${name%.fastq}`
     echo $name
     python3 $metalign_py $file  $metalign_data --output output_mergedR12_${name}.tsv --keep_temp_files --temp_dir ./temp_${name}
     unset file name
   done
   ```

   - [results](https://drive.google.com/open?id=1HFTW19Hl5OCHsoXQvd7PmbOOr6fNwFdc)

3. Validate the top novel species (2292357) by alignment

   - Code

   ```bash
   ### pull the corresponding genome out (code from Dr. David)
   cd /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_metalign_application/results/20200215_application1_variation/merged_R1R2_output/combined/temp_combined
   species=2292357
   rm pulled_species_${species}.fa 2>/dev/null 
   for accn in `grep -w $species subset_db_info.txt  | cut -f 1` 
   do
   	echo $accn
   	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < cmashed_db.fna | tail -n +2 | grep -A1 ${accn} >> pulled_species_${species}.fa
   	unset accn
   done
   
   ### build BWA and align the data
   original_fq="/gpfs/group/dmk333/default/shaopeng/projects/202002_metalign_application/results/20200215_application1_variation/merged_R1R2_output/original/merged_R12_Tomato_Pig_14_Day_14_S68_L002_original_001_original.fastq"
   reseq_fq="/gpfs/group/dmk333/default/shaopeng/projects/202002_metalign_application/results/20200215_application1_variation/merged_R1R2_output/reseq/merged_R12_Tomato_Pig_14_Day_14_S85_L001_reseq_001_resequence.fastq"
   combined_fq="/gpfs/group/dmk333/default/shaopeng/projects/202002_metalign_application/results/20200215_application1_variation/merged_R1R2_output/combined/merged_R12_Tomato_Pig_14_Day_14_combined.fastq"
   
   cd /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_metalign_application/results/20200224_validate_some_species
   # create index 
   fasta_file=pulled_species_2292357.fa
   echo "Creating BWA reference from fasta file"
   mkdir bwa_ref_${fasta_file}
   mv ${fasta_file} ./bwa_ref_${fasta_file}
   cd ./bwa_ref_${fasta_file}
   ml bwa
   bwa index ${fasta_file} \
     && echo "creating index successful!"
   # align data
   bwa_index=`find $PWD -maxdepth 2 -mindepth 2 -name "*.fa" ! -type d`
   time_tag=`date +"%H%M"`
   mkdir align_result_${time_tag}
   cd align_result_${time_tag}
   bwa mem ${bwa_index} ${original_fq} > aln_original.sam
   bwa mem ${bwa_index} ${reseq_fq} > aln_reseq.sam
   bwa mem ${bwa_index} ${combined_fq} > aln_combined.sam
   # collect mapping status
   for file in `ls *sam`
   do
     echo ${file}
     name=`echo ${file%.bam}`
     name=`echo ${name#aln_}`
     samtools flagstat ${file} > mapping_status_${name}.txt
     unset name file
   done
   ```

   - results:
     - /storage/home/sml6467/shaopeng_Koslicki_group/projects/202002_metalign_application/results/20200224_validate_some_species/align_result_1517
     - mapped ratio in C/O/R is 3.15% / 3.53% / 2.72%
     - Should be true signal, but only detected at combined data (double the depth)



### Run all data: <a name="runall"></a>

### Summary: <a name="summary"></a>



