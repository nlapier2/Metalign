# Metalign application 1 final report

Name: Shaopeng Liu

Email: sml6467@psu.edu

Last update: 2020.2.24

---

### Description:

1. Question: is it better to merge 2 replicates together or run them separately?

2. Sub-questions and conclusion:

   - Are 2 replicates' results consistent? **Yes**
   - Is result of merged data consistent with the results of 2 replicates? **Yes**
   - Which one is better in regard of finding the presence of species? **Merged data**

3. Summary:

   Metalign has low specificity (i.e. some false negative) for species with small percentage. One possible reason is that the low depth of species can only capture partial of its original genome, contributing to a low overlap to the known reference. So increasing sequencing depth can improve the capture of low-presence species. 

---

### Analysis:

O for original data; R for resequence data; C for combined of O and R.

1. Use metalign to process O/R/C R1 (read1) data, check the relationship
2. Use metalign to process O/R/C of R1+R2 (merge read1 and read2 directly) data, check the relationship. In this case, the depth is doubled then 1.
3. For the extra species identified in C, use alignment method (which can truely reflect presence) to validate if they are true signal.

----

### Detailed results:

1. for R1 file of O/R/C:
   - [original vs reseq](https://drive.google.com/file/d/1DUzZwB4OxxFc2S-g0H4vXP6EXBrRfMVb/view)
   - [combined vs original+reseq](https://drive.google.com/file/d/1OPzYJWfABIjTKW0zB0ckLJ-ulRhkbc7b/view)
   - Conclusion: 
     - Captured species: C=47, O=21, R=21. Only 1 species identified by O/R is not captured by C. (see part3 for more details)
     - the reseq and original data are consistent, although there are 2 non-overlap species in each one
     - the combined data are still consistent with original and reseq, while it identifies many **new species**, which lead to the deviation from y=x line.
2. for merged R1R2 file of O/R/C (note this time the depth is doubled compared with 1)
   - [original vs reseq](https://drive.google.com/open?id=1AxdgI5MSoXop3b3KPT4FwhuK_91EOe0z)
   - [combined vs original+reseq](https://drive.google.com/open?id=1VJYlewNxfVefp07NmlDC916UEqSrrU08)
   - Conclusion:
     - Captured species: C=76, O=36, R=33. All O and R are captured by C. (see part3 for more details)
     - similar to 1, the reseq and original are consistent, and the combined data is also consistent with the rest two.
     - When we double the depth (R1+R2), we found 62% more species in C, 71% more in O, and 57% more in R. Indicating that higher depth can capture more species.
3. [summary table for 1 and 2](https://drive.google.com/drive/folders/1HFTW19Hl5OCHsoXQvd7PmbOOr6fNwFdc)
4. use alignment method (BWA MEM) to validate the new identified species
   - Picked species ID: "2292357", this is the highest species in C that are not in O/R. Percentage is around 1.75%.
   - Align the merged R1R2 reads of C/O/R back to this species to check the presence of reads in the original data.
   - Results: C=3.15%, O=3.53%, R=2.72%
   - Conclusion:
     - The new species caputured in C is **true signal**. But it doesn't show up in O/R.
     - The percentage of reads belong to this species are similar in all C/O/R, so the reason why it can be captured by C is the depth.
5. All the scripts to reproduce the results above is recorded in file "application1_record.md".



