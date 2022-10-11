
## Algorithm idea

### Step (0) 
Take three input FASTAs (e.g. HG002, HG004, HG005)
(0B) Make sure FASTA seq IDs are both unique and have the sample name. Proposal: "SampleName_i_existingFastaSeqID" whereby i is the increment
(0C) Align these three FASTAs to the T2T reference, and create a BAM file

### Step (1) 
Calculate coverage in BAM file using a rolling window of N bp whereby N=100.
The output of this coverage file should be a tab-delimited BED-like file for each region whereby the coverage is >0
i.e.
```
START    END   COVERAGE
```

### Step (2) 
Download a "gene regions" BED file with [tss, end] coordinates for each gene


### Step (3) 
Calculate Bedtools intersect between these two BED files. Output the intersect BED files. These windows should have a unique name.

### Step (4)
For each [tss, end] window, calculate the bins starts/ends  (as @Yutong Qiu describes)

### Step (5)
For each window, extract the "reads" i.e. the isoforms which aligned in that window. We can output this into a separate FASTA file for each window. We'll need to name the output FASTA with the unique name given in Step (3)

### Step (6)
We now have sequences from three samples (HG002, HG004, HG005) for each window.

We will try doing intersections of sets first, based on the number of combinations.
e.g. "sample1 vs sample2", "sample2 vs sample 3", "sample1 vs sample3"

Ideas how to do this in python here

If there are intersections, we now have candidate unique isoforms on a sample basis.

Next steps:
Pairwise alignment -> variants.

