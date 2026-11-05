# Practical 4: Measuring differential expression of Transposable Elements with `TEtranscripts`

**Tutorial Author:** Clément Goubert (`goubert.clement@gmail.com`)

**Last Update**: November 4, 2026 (`TEtranscripts` 2.2.4, `STAR` 2.7.10b)

<!-- TODO(server, 5 Nov): nobody has run the commands below with /data/practical4_tetranscripts and the shared STAR index yet. Run them on the server, then update the date, the times, and the example outputs marked with TODO. -->

## Description

In this practical session, we will perform a differential expression analysis (statistical measurement of change in mRNA levels between conditions) of *Drosophila melanogaster*'s TEs after knockdown (KD) of Piwi,
the main gene involved in the Piwi-associated silencing complex (piRISC) which silences TEs in the Drosophila ovary.

![piRISC](https://i.imgur.com/brFjky5.png)
> Figure 1. The pi-RISC complex in *D. melanogaster*, adapted from [Onishi et al., 2021](https://doi.org/10.15252/embr.202153062)

This tutorial is adapted from the demonstration of `TEtranscripts` by [Jin and Hammell (2018)](https://link.springer.com/book/10.1007/978-1-4939-7710-9). Since then, `DESeq2`, the program that performs the differential expression analysis, requires replicates for each condition. Thus, we have "upgraded" the analysis to use a new dataset from [Ilyin et al., 2021](https://www.mdpi.com/1422-0067/22/24/13430), as well as more recent reference genome (dm6) and TE annotations.

In their experiments, Ilyin et al., investigate the interaction between Piwi and HP1a for the repression of TEs in the ovaries of *D. melanogaster*. For this tutorial, we will only focus on the Piwi knock-down (Piwi-KD) vs control datasets, in order to see what happens when this key protein is not produced.

The main challenge with TE RNA-seq is the abundance of ambiguous reads which can map equally to different TE copies. This challenge is illustrated in the following figure:

![Mappability](https://i.imgur.com/mSfw57Y.png)
> Figure 2. Origin of multi-mapping reads in TE-RNAseq. From [Lanciano and Cristofari (2020)](https://www.nature.com/articles/s41576-020-0251-y).

Different strategies exist to map and quantify TE expression, and I encourage those who would like to know more about it to read [this 2020 review by Lanciano and Cristofari](https://www.nature.com/articles/s41576-020-0251-y). We will discuss the solution used by the software `TEtranscripts` which uses a reference-genome approach (vs. reference transcriptome or hybrid approach). The key idea is to first identify multi-mapping reads, and then attempt to assign them to a unique locus. If the read can't be assigned unambiguously, then its weight (w = 1), will be equally distributed to all competing (n) loci (w = 1/n). Then, an Expectation-Maximization algorithm, will iteratively use information from uniquely mapped reads to find the most likely loci (or best set of loci) of origin. In the end, `TEtranscripts` will sum the weights (normalized counts) of uniquely and multi-mapping reads to the most suitable TE family (Family-level quantification). This family-level approach has been proven to be reliable, with high reproducibility. Note that other methods such as [SQuIRE](https://academic.oup.com/nar/article/47/5/e27/5280934) or [TElescope](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006453) have been developed to estimate TE abundances via EM-algorithms at the locus level. Another notable method is [Cello-seq](https://www.nature.com/articles/s41587-021-01093-1), which allows TE DE analysis at the copy level AND in single-cell data!

![EM](https://i.imgur.com/Ih9H46n.png)
> Figure 3. Simplified representation of algorithms using read-mapping to a reference genome. The red boxes mark the parts specific to `TEtranscripts`. Adapted from [Lanciano and Cristofari (2020)](https://www.nature.com/articles/s41576-020-0251-y).

### Objectives

1. Align RNA-seq reads to the *D. melanogaster* reference genome with `STAR` aligner
2. Estimation of relative transcript abundance with `TEtranscripts`
3. Statistical Differential Expression (DE) analysis with `TEtranscripts`

## Protocol

### Working directory set-up

`STAR`, `TEtranscripts` and `R` with `DESeq2` are in the conda environment `tetranscripts`. Activate it first:

```
conda activate tetranscripts
```

We set up a working directory that will host our results

```
mkdir practical4
cd practical4
```
Next, we create a symbolic link to the input data. This allows us to work as if the data were located in the same working directory (like in real life!).

```
ln -s /data/practical4_tetranscripts/input_data/ input_data
```

### Input data

Let's first look at the input data, located in the shared folder of the course server

```
tree input_data/
```

<!-- TODO(server, 5 Nov): paste the real output of tree here. -->

```
input_data/
├── Control_rep_1.fastq.gz      # --> control RNA-seq reads, replicate 1
├── Control_rep_2.fastq.gz      # --> control RNA-seq reads, replicate 2
├── dm6.ensGene.gtf             # --> dm6 genes annotation
├── dm6.fa                      # --> dm6 reference genome
├── dm6_rmsk_TE.gtf             # --> dm6 TE annotation
├── Piwi_GKD_rep_1.fastq.gz     # --> Piwi KD RNA-seq reads, replicate 1
└── Piwi_GKD_rep_2.fastq.gz     # --> Piwi KD RNA-seq reads, replicate 2
```
>Note: for future use, check details about the `gtf` files required for `TEtranscripts` [here](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts), and a way to convert `.gff` to `.gtf` [here](https://www.biostars.org/p/45791/)

### Step 1: Align RNA-seq to the reference genome

#### 1.1 The STAR index for the dm6 reference genome

`STAR` first needs an index of the reference genome. With `--sjdbGTFfile`, the index also stores the splice junctions of the gene annotation. Building the index took about 7 minutes and 5.5 GB of memory in our tests. With 20 students building it at once, the course server would run out of memory, so we built one index for everyone with this command (**do not run it**):

```
STAR --runMode genomeGenerate --genomeDir Index --genomeFastaFiles input_data/dm6.fa --sjdbGTFfile input_data/dm6.ensGene.gtf --sjdbOverhang 100 --runThreadN 2 --genomeSAindexNbases 12
```

Link the shared index into your working directory instead:

```
ln -s /data/practical4_tetranscripts/Index/ Index
```

#### 1.2 align the reads

With the index created, we can now **map the RNA-seq reads for the control and KD conditions**. Because `TEtranscripts` relies on differences between **uniquely mapping** and **multi-mapping reads**, we will use specific parameters of the `STAR` aligner to allow it (`--winAnchorMultimapNmax` and `--outFilterMultimapNmax`). The index already holds the splice junctions of the gene annotation, so the alignment commands do not need the `gtf` file.

#### 1.2.1 Treatment (2 replicates)
```
mkdir STARalignments
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Piwi_GKD_rep_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiKD1_
```
>note: each alignment should take about 1 minute <!-- TODO(server, 5 Nov): time and memory per alignment without --sjdbGTFfile -->

```
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Piwi_GKD_rep_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiKD2_
```

#### 1.2.2 Control (2 replicates)
```
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Control_rep_1.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiCT1_
```
```
STAR --runThreadN 2 --genomeDir Index --readFilesIn input_data/Control_rep_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --winAnchorMultimapNmax 200 --outFilterMultimapNmax 100 --outFileNamePrefix STARalignments/piwiCT2_
```

### Step 2: Run TEtranscripts to perform counts and differential expression analysis

```
TEtranscripts --format BAM --stranded reverse -t STARalignments/piwiKD1_Aligned.out.bam STARalignments/piwiKD2_Aligned.out.bam  -c STARalignments/piwiCT1_Aligned.out.bam STARalignments/piwiCT2_Aligned.out.bam --GTF input_data/dm6.ensGene.gtf --TE input_data/dm6_rmsk_TE.gtf --mode multi --project piwiKD_vs_control --minread 1 -i 10 --padj 0.05
```

> The message: `Could not retrieve index file for 'STARalignments/piwiKD1_Aligned.out.bam'` will appear. It is not an error, `TEtranscripts` will actually generate the index for each alignment automatically.
> `TEtranscripts` should take about 6 minutes to run. <!-- TODO(server, 5 Nov): check the message and the time -->

While `TEtranscripts` is running, you can look into the detailed workflow. You will see that the program alternated between gathering reads from each alignments and performing the Expectation/Maximization step attributing reads to TE loci.

![TEtranscripts](https://i.imgur.com/BKufoQh.gif)
> Figure 4. The `TEtranscripts` flow chart, from Figure 1 of [Jin et al., 2015](https://doi.org/10.1093/bioinformatics/btv422).


### Results

After a successful run, 4 files will be generated:

<!-- TODO(server, 5 Nov): the four example outputs below come from the 2023 handout. Replace them with the output of the 5 Nov run; the first line of the cntTable is a header (gene/TE, then one column per BAM file). -->

- `piwiKD_vs_control.cntTable`: this is the read counts out of `TEtranscripts`, for both genes and TEs. The first column indicate gene or TE family name, and each next column correspond to the read count for each condition.

    ```
    head piwiKD_vs_control.cntTable
    ```

    ```
    "FBgn0000003"   9067    0   15412   5975
    "FBgn0000008"   32  28  47  59
    "FBgn0000014"   80  75  25  19
    "FBgn0000015"   8   5   1   2
    "FBgn0000017"   158 207 115 153
    "FBgn0000018"   27  21  29  27
    "FBgn0000022"   1   0   0   0
    "FBgn0000024"   7   8   2   2
    "FBgn0000028"   1   0   0   0
    ```
    >*D. melanogaster* genes all have the `FBgn` prefix
    ```
    tail piwiKD_vs_control.cntTable
    ```

    ```
    TV1I:Gypsy:LTR  5   8   7   5
    Tn1000:ARTEFACT:ARTEFACT    0   0   0   0
    Transib-N1_DM:CMC-Transib:DNA   2   2   1   3
    Transib5:CMC-Transib:DNA    10  9   9   13
    UHU:TcMar-Tc1:DNA   0   0   0   0
    Ulysses_I-int:Gypsy:LTR 2   5   1   0
    XDMR:XDMR:Unknown   0   2   2   4
    XDMR_DM:XDMR_DM:Unknown 3   1   0   0
    ZAM_I-int:Gypsy:LTR 8   9   4   6
    ZAM_LTR:Gypsy:LTR   0   2   6   6
    ```

- `piwiKD_vs_control_gene_TE_analysis.txt`: this is the main result table from the differential expression analysis with `DESeq2`.

    ```
    head piwiKD_vs_control_gene_TE_analysis.txt
    ```

    ```
    baseMean    log2FoldChange  lfcSE   stat    pvalue  padj
    FBgn0000003 7496.32149464016    -1.31624683855759   3.91511703871456    -0.336196038468815  0.736723040072146   0.90321224259372
    FBgn0000008 41.2452052205847    -0.846404602618286  0.468447081016323   -1.80683077538227   0.0707886686420089  0.243117726171184
    FBgn0000014 49.1262294684561    1.79266376608722    0.453418090720179   3.95366616986955    7.69627701796716e-05    0.00110998817889254
    FBgn0000015 3.92336261262419    2.08258292784614    1.4980808020054 1.39016728941343    0.164478083890275   0.421438277712362
    FBgn0000017 157.361890606094    0.429714476761589   0.309100895585062   1.39020780236897    0.16446578486819    0.421438277712362
    FBgn0000018 25.7681516759187    -0.250195946949605  0.565873845705727   -0.442140856744446  0.65838728075958    0.865750510693057
    FBgn0000024 4.70352527222641    1.88724067149184    1.34287753608743    1.40537064682044    0.159911090558418   0.415032094465402
    FBgn0000032 50.7878534373474    -0.532402717165373  0.432720270103867   -1.23036232399647   0.2185614553133 0.503891826335278
    FBgn0000042 968.914785124894    0.248341515484232   0.198409438947303   1.25166179997208    0.210693125186219   0.493176903629541
    ```

    ```
    tail piwiKD_vs_control_gene_TE_analysis.txt
    ```

    ```
    TRANSPAC_I-int:Gypsy:LTR    116.414503755064    2.06051339821219    0.331614247696445   6.21358525010769    5.17891139816696e-10    3.41974291063766e-08
    TRANSPAC_LTR:Gypsy:LTR  12.108452499695 1.76594902622038    0.853463486837924   2.06915591991312    0.0385314576487211  0.159105911832808
    TV1I:Gypsy:LTR  6.23043397426447    0.102040532049877   1.10622965622402    0.0922417252835001  0.926505987354992   0.975080697666644
    Transib-N1_DM:CMC-Transib:DNA   1.98525372581525    -0.0227181921579707 1.91002621820881    -0.0118941781748291 0.990510042631103   NA
    Transib5:CMC-Transib:DNA    10.1709658430894    -0.236336499485258  0.871958330375934   -0.271041047779616  0.786359461879388   0.928960553907422
    Ulysses_I-int:Gypsy:LTR 2.00205880394992    2.80416164448381    2.28677817067709    1.2262499618201 0.220104625120975   NA
    XDMR:XDMR:Unknown   2.01165924819134    -1.57586245854514   2.07183136671284    -0.760613283428274  0.446888083329893   NA
    XDMR_DM:XDMR_DM:Unknown 0.966253597048695   3.39638119897609    2.9987731534211 1.13259023781155    0.257386370097645   NA
    ZAM_I-int:Gypsy:LTR 6.69641232303531    0.745381479347682   1.07298376376424    0.694681042267348   0.48725521121973    0.764045968153913
    ZAM_LTR:Gypsy:LTR   3.5121624278533 -2.5835928062637    1.6998712045284 -1.51987562315375   0.128542238135249   0.362461822414429
    ```
    >REMINDER: with `TEtranscripts`, the results are given at the TE family level

- `piwiKD_vs_control_sigdiff_gene_TE.txt` is the same table, filtered for significant genes or TE families only (according to the threshold set in the `TEtranscripts` command)

- `piwiKD_vs_control_DESeq2.R` is a `R` script that can be used to replicate the DE analysis, changing parameters, etc... Without re-running the read counting step of `TEtranscripts`


#### Plot the results as a volcano plot!

Before analyzing the results, create a symlink to the actual results, made out of the complete read-sets.

```
ln -s /data/practical4_tetranscripts/results/ results
```

A few guidelines:
- Volcano plot represent the relationship between log2-fold-change (log2FC), i.e. the intensity of the change, and the significance of the change (P-value)
- We will use the adjusted P-value (Padj < 0.05) for our P-value significance threshold
    - It is more handy to represent the adjusted P-value as -log10(Padj), such as the more significant, the higher the value.
- We will use |log2FC| > 2 (a change of more than 4 times, up or down) for our fold-change significance. This threshold is empirically supported by qPCR validation.

Use the table `results/piwiKD_vs_control_gene_TE_analysis.txt` and the tool of your choice to generate volcano plots. Below is an example made with base `R`; the code is at the end of this section. Excel is more than OK!

To transfer files from the server to your computer:

1. Open a new terminal tab or window
2. From this new tab, copy the result table to the folder you are in on your computer:

    ```
    scp -i <KEY>.pem <USER>@<SERVER_IP>:practical4/results/piwiKD_vs_control_gene_TE_analysis.txt .
    ```

    Replace `<KEY>.pem` with your key file, `<USER>` with your user name on the server, and `<SERVER_IP>` with the server address you receive on the day of the class. You can then open the table from your file explorer.
3. Alternatively, you can use FileZilla to transfer the output table.

![volcano](figures/volcano_piwiKD.png)
> Figure 5. Volcano plots of the Piwi-KD vs control comparison. Left: all genes and TE families. Right: TE families only; in red, those with Padj < 0.05 and |log2FC| > 2. Made with `VolcanoPiwi.R` from the results in `results/`.

As expected, upon KD of Piwi, several TE families became transcriptionally active in the ovaries!

The `R` code for these plots is in `VolcanoPiwi.R`. To run it on the server, from your `practical4` directory:

```
cp /data/practical4_tetranscripts/VolcanoPiwi.R .
Rscript VolcanoPiwi.R
```

`Rscript` writes the two plots to `Rplots.pdf`. Copy it to your computer with `scp`, as above, using `practical4/Rplots.pdf` as the file to copy.

🎉 **Success!!!** Now you are an expert in TE differential expression analysis!
