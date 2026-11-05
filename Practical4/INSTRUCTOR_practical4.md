# Practical 4 (TEtranscripts): instructor run sheet

The commands to type live, in order. They are the commands of the student handout, `Practical4_TEtranscripts.md`.

<!-- TODO(server, 5 Nov): the times below come from the test account tetest on the earlier m4.xlarge (4 CPUs, 15.6 GB), with STAR building its own index and the alignments run with --sjdbGTFfile. Replace them with times measured on the m4.4xlarge with the shared index. -->

We measured the times as the test account `tetest`, one run at a time. We have not measured a run with 20 students at once; expect it to be slower.

## Before the class

- Get the server's address for the day (it changes every session).
- The data, the results made from the complete read sets, `VolcanoPiwi.R` and the shared STAR index are in `/data/practical4_tetranscripts`, read-only for students.
- Students do not build the STAR index. Building it took 395 s and 5.5 GB of memory; 20 builds at once would need about 110 GB, and the server has 62 GB.

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

## 2. STAR alignments (about 1 minute each, 4 minutes in all)

```
mkdir STARalignments
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Piwi_GKD_rep_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiKD1_
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Piwi_GKD_rep_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiKD2_
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Control_rep_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiCT1_
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Control_rep_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiCT2_
```

Each alignment took about 50 s and 3.7 GB of memory with the 2023 commands, which passed `--sjdbGTFfile` and used the account's own index. At 3.7 GB each, 20 students aligning at the same moment would need about 74 GB. <!-- TODO(server, 5 Nov): peak memory without --sjdbGTFfile, with the shared index -->

To see how many reads mapped once or to several places:

```
grep -E 'Uniquely mapped reads %|% of reads mapped to multiple loci' STARalignments/*_Log.final.out
```

## 3. TEtranscripts (about 6 minutes)

```
TEtranscripts --format BAM --stranded reverse -t STARalignments/piwiKD1_Aligned.out.bam STARalignments/piwiKD2_Aligned.out.bam  -c STARalignments/piwiCT1_Aligned.out.bam STARalignments/piwiCT2_Aligned.out.bam --GTF input_data/dm6.ensGene.gtf --TE input_data/dm6_rmsk_TE.gtf --mode multi --project piwiKD_vs_control --minread 1 -i 10 --padj 0.05
```

It took 345 s and 0.76 GB of memory. While it runs, go through Figure 4 of the handout (the `TEtranscripts` flow chart).

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

`Rscript` writes `Rplots.pdf`, the plot in Figure 5 of the handout. In the copy of `results/` in Google Drive, 11 TE families have Padj < 0.05 and |log2FC| > 2: 10 up in Piwi-KD and 1 down (`Gypsy3_LTR:Gypsy:LTR`).

## 6. Clean up (seconds)

The four BAM files should take about 400 MB per student (estimated from the 2023 handout, not measured). <!-- TODO(server, 5 Nov): check the size --> To remove them:

```
rm STARalignments/*_Aligned.out.bam
```
