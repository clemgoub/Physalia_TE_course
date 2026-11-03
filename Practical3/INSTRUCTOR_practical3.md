# Practical 3 (GraffiTE): instructor run sheet

The commands to type live as `ubuntu`, in order, with how long each took on 3 Nov 2026. They are the commands of the student handout, `Practical3_GraffiTE_v1.1.md`.

All times come from runs on the course server (m4.4xlarge, 16 CPUs, 62 GB) with nobody else running GraffiTE. The full run took 17 min 42 s as `ubuntu` and 17 min 42 s as the student test account `tetest`. We have not measured a run with 20 students at once; expect it to be slower.

## Before the class

- Get the server's address for the day (it changes every session).
- A finished run of the same command is in `/home/ubuntu/practical3_precomputed`, and a read-only copy of its `tutorial/` folder is in `/data/practical3_graffite/precomputed/tutorial`, for everyone.
- The results of a run without `--human` are in `/data/practical3_graffite/run_without_human`. Students look at them while their own run goes.

## 1. Set up (under 1 minute)

```bash
ssh -i <KEY>.pem ubuntu@<SERVER_IP>
```

Replace `<KEY>.pem` with the key file and `<SERVER_IP>` with the server address you receive on the day of the class.

```bash
mkdir practical3
cd practical3
ln -s /data/practical3_graffite/input/ input
tree input/
cat input/HG002.assemblies.csv
cat input/HG002.reads.csv
```

Each command returns at once. GraffiTE needs no conda environment, so `ubuntu`'s own conda setup does not matter here.

## 2. Start GraffiTE (17 min 42 s)

```bash
nextflow run /opt/graffite/main.nf \
  --assemblies input/HG002.assemblies.csv \
  --reference input/chm13v2.0.chr22.fa \
  --TE_library input/Homo_sapiens_DFAM3.9_040926.fa \
  --genotype_with input/HG002.reads.csv \
  --graph_method giraffe --human \
  --out tutorial \
  -qs 1 --cores 2 --stSort_m 1G --stSort_t 2 --make_graph_memory 4G
```

What the run was doing, in minutes after you press Enter (from the `ubuntu` run's Nextflow trace):

| From | To | Step | Process on screen |
|---|---|---|---|
| 0:00 | 9:20 | 1 | `map_asm`: aligns the maternal assembly, then the paternal one (4 min 38 s each) |
| 9:20 | 9:45 | 1 | `svim_asm`, `truvari_merge`: calls and merges the SVs |
| 9:45 | 11:25 | 2 | `repeatmask_VCF`: RepeatMasker on the SV sequences |
| 11:25 | 12:55 | 3 | `tsd_search`: TSD search, 10 batches of about 10 s |
| 12:55 | 15:50 | 4 | `concat_repeatmask`, then `make_graph`: builds the graph (2 min 50 s) |
| 15:50 | 16:35 | 4 | `hervk_annotate`, then `graph_align_reads`: aligns the reads with Giraffe |
| 16:35 | 17:40 | 4 | `vg_call`, then `merge_VCFs`, `genotyping_audit`, `hervk_reconcile` |

To skip to the results at any point, stop the run with `Ctrl + C` and have everyone type, in their `practical3` directory (each command returns at once):

```bash
mv tutorial tutorial_interrupted
ln -s /data/practical3_graffite/precomputed/tutorial tutorial
```

If `mv` says that `tutorial` does not exist, the `ln` command alone is enough. We checked on 3 Nov 2026, as the student test account, that every command of section 4 gives the same output on the shared copy.

At the end, Nextflow prints `Succeeded` and the run time. `tutorial/` then holds 838 MB, and `work/` 928 MB.

## 3. While it runs: results for any species (each command under 1 second)

```bash
ls /data/practical3_graffite/run_without_human/*
```

```bash
grep '##bcftools_viewCommand' /data/practical3_graffite/run_without_human/3_TSD_search/pangenome.trusted.vcf | tail -1
```

```bash
awk -F'\t' 'NR == 1 || $22 != "NA"' /data/practical3_graffite/run_without_human/4_Genotyping/GraffiTE.merged.genotypes.presence-absence_trusted.tsv | \
  cut -f 1,2,4,6,9,10,20,21,22 | head -8 | column -t
```

```bash
awk -F'\t' 'NR > 1 {print $10, $22}' /data/practical3_graffite/run_without_human/4_Genotyping/GraffiTE.merged.genotypes.presence-absence_trusted.tsv | sort | uniq -c
```

```bash
awk -F'\t' 'NR == 1 || $10 ~ /^(LINE|SINE|LTR|DNA|RC|Retroposon)/' \
  /data/practical3_graffite/run_without_human/4_Genotyping/GraffiTE.merged.genotypes.presence-absence_trusted.tsv > trusted_TEs.tsv
```

The fourth command shows the trusted subset is mostly satellites (336 of 368 records); the last one keeps the TE classes only.

## 4. After the run: outputs (each command under 1 second)

```bash
tree -L 2 tutorial/
grep -v '#' tutorial/1_SV_search/SVs.vcf | cut -f 3 | grep -o 'svim_asm\.[A-Z]*' | sort | uniq -c
grep -v '#' tutorial/3_TSD_search/pangenome.vcf | cut -f 3 | grep -o 'svim_asm\.[A-Z]*' | sort | uniq -c
grep -v '#' tutorial/1_SV_search/SVs.vcf | grep -w 'HG002.mat.svim_asm.DEL.71_275'
grep -v '#' tutorial/3_TSD_search/pangenome.vcf | grep -w 'HG002.mat.svim_asm.DEL.71_275' | cut -f 1-3,6-11
```

```bash
zcat tutorial/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz | grep '^#CHROM'
zcat tutorial/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz | grep -v '^#' | cut -f 10 | cut -d: -f 1 | sort | uniq -c
zcat tutorial/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz | grep -w 'HG002.pat.svim_asm.DEL.140_387' | cut -f 1-3,6-10
zcat tutorial/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz | grep -w 'HG002.mat.svim_asm.INS.182_786' | cut -f 1-3,6-10
grep -v '#' tutorial/3_TSD_search/pangenome.vcf | grep -w 'HG002.mat.svim_asm.DEL.83_330' | cut -f 1-3,8
```

The human subset:

```bash
sed -n '/^records in/,/^dropped/p' tutorial/3_TSD_search/human_filter_summary.txt
cut -f 1,2,4,6,9,10,22,23 tutorial/3_TSD_search/pangenome.presence-absence_human.tsv | head -6 | column -t
```

The summary per superfamily, with SVA VNTR-only records apart from SVA insertions:

```bash
zcat tutorial/4_Genotyping/GraffiTE.merged.genotypes.human.vcf.gz | grep -v '^#' | \
  awk -F'\t' '
  { match($8, /matching_classes=[^;]*/); sf = substr($8, RSTART + 17, RLENGTH - 17)
    if ($8 ~ /VNTR_only/) sf = "SVA_VNTR_only"
    match($8, /SVLEN=[^;]*/); len = substr($8, RSTART + 6, RLENGTH - 6) + 0
    split($10, f, ":"); gt = f[1]
    if (gt ~ /\./) col = "no_call"
    else { alt = gsub(/1/, "1", gt); col = (len > 0 ? alt : 2 - alt) " copies" }
    n[sf]++; t[sf, col]++ }
  END { print "superfamily", "records", "TE_on_0_copies", "TE_on_1_copy", "TE_on_2_copies", "no_call"
        for (sf in n) print sf, n[sf], t[sf, "0 copies"] + 0, t[sf, "1 copies"] + 0, t[sf, "2 copies"] + 0, t[sf, "no_call"] + 0 }' | \
  column -t
```

Expected result:

```
superfamily     records  TE_on_0_copies  TE_on_1_copy  TE_on_2_copies  no_call
LINE/L1         4        0               0             0               4
Retroposon/SVA  2        1               0             0               1
SVA_VNTR_only   1        0               0             0               1
SINE/Alu        21       9               4             1               7
```

`awk` prints the rows in no fixed order, so they may come out in a different order on screen.

## 5. Clean up (seconds)

```bash
rm -rf work/
```

On the shared results there is no `work/` to remove.
