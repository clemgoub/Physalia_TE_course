# Practical 4 (TEtranscripts): instructor run sheet

The commands to type live, in order. They are the commands of the student handout, `Practical4_TEtranscripts.md`.

We measured the times on 5 November 2026 as the test account `tetest` on the course server (16 CPUs, 62 GB of memory), one run at a time. The steps we timed add up to about 6 minutes. We have not measured a run with 20 students at once; expect it to be slower.

## Before the class

- Get the server's address for the day (it changes every session).
- The data, the results made from the complete read sets, `VolcanoPiwi.R` and the shared STAR index are in `/data/practical4_tetranscripts`, read-only for students.
- Students do not build the STAR index. With the handout's command (2 threads), building it took 395 s and 5.5 GB of memory on the m4.xlarge we tested on earlier (4 CPUs, 15.6 GB). 20 builds at once would need about 110 GB, and the server has 62 GB. We built the shared index in 92 s with 8 threads; it used 5.2 GB.

## 1. Set up (under 1 minute)

```
ssh -i <KEY>.pem ubuntu@<SERVER_IP>
```

Replace `<KEY>.pem` with the key file and `<SERVER_IP>` with the server address you receive on the day of the class.

```
conda activate tetranscripts
mkdir practical4
cd practical4
ln -s /data/practical4_tetranscripts/input_data/ input_data
tree input_data/
ln -s /data/practical4_tetranscripts/Index/ Index
```

Show the index command from step 1.1 of the handout, but do not run it.

## 2. STAR alignments (about 30 seconds each, 2 minutes in all)

```
mkdir STARalignments
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Piwi_GKD_rep_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiKD1_
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Piwi_GKD_rep_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiKD2_
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Control_rep_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiCT1_
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Control_rep_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiCT2_
```

Each alignment took 27 to 29 s and 2.3 GB of memory. At 2.3 GB each, 20 students aligning at the same moment would need about 45 GB (not measured).

The 2023 commands also passed `--sjdbGTFfile` and `--sjdbOverhang` to each alignment. The shared index already holds the junctions from the GTF. We aligned `Piwi_GKD_rep_1` with and without these two options and got the same SAM records and the same `Log.final.out` statistics; with them, the alignment took 37 s and 3.7 GB instead of 29 s and 2.3 GB.

To see how many reads mapped once or to several places:

```
grep -E 'Uniquely mapped reads %|% of reads mapped to multiple loci' STARalignments/*_Log.final.out
```

On 5 November, the percentages were:

| Sample | Mapped once | Mapped to several places |
|---|---|---|
| piwiKD1 | 76.28% | 22.75% |
| piwiKD2 | 71.53% | 27.33% |
| piwiCT1 | 78.26% | 20.94% |
| piwiCT2 | 80.16% | 18.73% |

## 3. TEtranscripts (about 4 minutes)

```
TEtranscripts --format BAM --stranded reverse -t STARalignments/piwiKD1_Aligned.out.bam STARalignments/piwiKD2_Aligned.out.bam  -c STARalignments/piwiCT1_Aligned.out.bam STARalignments/piwiCT2_Aligned.out.bam --GTF input_data/dm6.ensGene.gtf --TE input_data/dm6_rmsk_TE.gtf --mode multi --project piwiKD_vs_control --minread 1 -i 10 --padj 0.05
```

It took 250 s and 0.76 GB of memory. It prints `Could not retrieve index file` once for each BAM file; this is not an error. While it runs, go through Figure 4 of the handout (the `TEtranscripts` flow chart).

## 4. Results (each command under 1 second)

```
head piwiKD_vs_control.cntTable
tail piwiKD_vs_control.cntTable
head piwiKD_vs_control_gene_TE_analysis.txt
tail piwiKD_vs_control_gene_TE_analysis.txt
ln -s /data/practical4_tetranscripts/results/ results
```

## 5. Volcano plot (seconds)

```
cp /data/practical4_tetranscripts/VolcanoPiwi.R .
Rscript VolcanoPiwi.R
```

`Rscript` writes `Rplots.pdf`, the plot in Figure 5 of the handout. The files in `results/` on the server are the same as the copy in Google Drive (same md5). In them, 11 TE families have Padj < 0.05 and |log2FC| > 2: 10 up in Piwi-KD and 1 down (`Gypsy3_LTR:Gypsy:LTR`).

## 6. Clean up (seconds)

The four BAM files take 321 MB per student, out of 325 MB for the whole `practical4` directory. To remove them:

```
rm STARalignments/*_Aligned.out.bam
```
