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
     - [ ] follow1: align 2 reps to the extra species in merged file (5%)
   - [ ] Compare the merged R1R2 of the 2 rep of "Pig_14_Day_14" data
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

   ```R
   # Check consistency in R
   
   ```

   

### Run all data: <a name="runall"></a>

### Summary: <a name="summary"></a>



