# purvis-dev2025-crispr-outcomes

Purvis Dev 2025 CRISPR editing outcomes analysis. It uses the `CrispRVariants` package to analyze CRISPR editing outcomes from paired-end sequencing data. Amplicon sequencing data is used to assess the efficiency of CRISPR editing at a specific target region in the mouse genome.

## Status: In Progress

Start Date: 2025-07-30

Hours spent:

## Project Notes

### Investigators

Joseph Brzezinski IV, Principal Investigator
<JOSEPH.BRZEZINSKI@CUANSCHUTZ.EDU>

Ian Purvis, Graduate Student
<IAN.PURVIS@CUANSCHUTZ.EDU>

### Background
This project aims to analyze CRISPR editing outcomes using paired-end sequencing data. The `CrispRVariants` package is utilized to assess the efficiency of CRISPR editing at a specific target region in the mouse genome. The analysis will focus on identifying and quantifying indels (insertions and deletions) resulting from CRISPR editing for the revision of the A robust cis-regulatory network ensures Otx2 expression during retinal development manuscript.

### Hypothesis

CRISPR editing will result in a higher frequency of indels at the target region compared to non-target regions.

### Samples

  P0 mouse retina samples electroporated with CRISPR NHEJ expressing plasmid reagents targeting the Otx2 gene or enhancer regions

| sampleid     | prefix              | treatment |
|--------------|---------------------|-----------|
| ControlDHS2  | ControlDHS2primers  | Control   |
| ControlDHS15 | ControlDHS15primers | Control   |
| TestDHS15    | TestDHS15primers    | Test      |


**guide RNA sequences used for CRISPR editing:**
  
DHS15- g1 = GCTGCCCCAGCCTTTCACAA
DHS15 - g2 = TTTTCTTTTTTATTTAACCG
DHS2 - g1 = CACCGACCCAGCTCTGGCAGATGGG
DHS2 - g2 = ATTTTGAACCCCCTACAGAT

These sequences do not include the PAM sequence.

### Raw Data

Analysis Copy:  
`https://github.com/MLKaufman/purvis-dev2025-crispr-outcomes/tree/main/raw_data`

Investigator Copy:  

---

## Analysis

- Process and align paired-end sequencing FASTQ files using `bwa` and `samtools`
- Use `CrispRVariants` to analyze CRISPR editing outcomes from the aligned BAM files
- Identify and quantify indels at the target region

### Notes

- mm39 genome build will be used for alignment

---

### Objectives

- Analyze CRISPR editing outcomes in mouse retina samples
- Identify and quantify indels at the target region

### Progress

- [x] Install `CrispRVariants` package and preliminary setup
- [x] Process paired-end sequencing FASTQ files
- [ ] Generate detailed mutagenesis analysis report
- [ ] Create coverage plots for each sample

### Deliverables
<https://github.com/MLKaufman/purvis-dev2025-crispr-outcomes/tree/main/analysis/outs/>

- [ ] Mutagenesis analysis report `01-CrispRVariants.html`
- [ ] files and/or coverage plots for each sample

### Apps

---

## Findings

### Summary

### Publication

### GEO Submission

---

## References / Support

https://www.bioconductor.org/packages/release/bioc/vignettes/CrispRVariants/inst/doc/user_guide.html

https://github.com/markrobinsonuzh/CrispRVariants

https://bioconductor.org/packages/release/bioc/html/CrispRVariants.html