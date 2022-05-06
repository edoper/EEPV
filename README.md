# EEPV
A framework to identify and visualize Exons Enriched with Pathogenic Variants

## Introduction
[Insertar abstract del paper aqui.]
## Contents
- **Content 1**: Asd.
- **Content 2**: Asd.
- **Content 3**: Asd.	

## Pipeline
The pipeline for EEPV starts with the downloaded GnomAD data (https://gnomad.broadinstitute.org/).
  ```
  gnomad/release/2.1.1/vcf/exomes/
  gnomad/release/2.1.1/vcf/genomes/
  ```
We require bcftools (https://samtools.github.io/bcftools/bcftools.html) and tabix (http://www.htslib.org/doc/tabix.html) for processing.

### Part 1: GnomAD
After downloading the data, we begin extracting the indels to make a consolidated vcf for further processing.
```
# Extracting indels from GnomAD exomes.
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
	echo "working with $K"
	bcftools view --types indels gnomad.exomes.r2.1.1.sites.$K.vcf.bgz >$K.indels.exome.vcf
	bgzip -c $K.indels.exome.vcf >$K.indels.exome.vcf.gz
	tabix -p vcf $K.indels.exome.vcf.gz
done
# Extracting indels from GnomAD genomes.
bcftools view --types indels gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz >all.indels.genome.vcf
bgzip -c all.indels.genome.vcf >all.indels.genome.vcf.gz
tabix -p vcf all.indels.genome.vcf.gz
```

Since the GnomAD exome data is divided by chromosomes, we merge it and create a final vcf file.
```
# Exome vcf merging
bcftools concat -O z -o all.indels.exome.vcf.gz 1.indels.exome.vcf.gz 2.indels.exome.vcf.gz 3.indels.exome.vcf.gz 4.indels.exome.vcf.gz 5.indels.exome.vcf.gz 6.indels.exome.vcf.gz 7.indels.exome.vcf.gz 8.indels.exome.vcf.gz 9.indels.exome.vcf.gz 10.indels.exome.vcf.gz 11.indels.exome.vcf.gz 12.indels.exome.vcf.gz 13.indels.exome.vcf.gz 14.indels.exome.vcf.gz 15.indels.exome.vcf.gz 16.indels.exome.vcf.gz 17.indels.exome.vcf.gz 18.indels.exome.vcf.gz 19.indels.exome.vcf.gz 20.indels.exome.vcf.gz 21.indels.exome.vcf.gz 22.indels.exome.vcf.gz X.indels.exome.vcf.gz Y.indels.exome.vcf.gz
tabix -p vcf all.indels.exome.vcf.gz
```

We filter the generated vcfs to keep only variants with the PASS label.
```
# PASS variant filtering.
bcftools view -O z -f 'PASS' all.indels.exome.vcf.gz >all.indels.exome.pass.vcf.gz
bcftools view -O z -f 'PASS' all.indels.genome.vcf.gz >all.indels.genome.pass.vcf.gz
tabix -p vcf all.indels.exome.pass.vcf.gz
tabix -p vcf all.indels.genome.pass.vcf.gz
```

Finally, we merge both genome and exome indels.
```
# Vcf merging.
bcftools merge -O z -o indels.vcf.gz all.indels.exome.pass.vcf.gz all.indels.genome.pass.vcf.gz
```









- **Download**: Download the repository and uncompress the db folder contents. Make sure you have installed the Perl modules “*Data::Dumper*” and “*List::MoreUtils*” and the R packages "*ggplot2*", "*readr*" and "*ggrepel*" before running the code. 
- **Command**: From terminal and inside the repository directory run:

  `$ perl Part-1-missense-aligner.pl`

- **Output**: In the db/ folder, two files will be produced per aligment, one of them accounting for clinvar-hgmd missense variant mapping over the alignment (family-name.clinvar-hgmd.binary) and the other accounting for gnomad missense variants mapping over the same alignment (family-name.gnomad.binary).
- **Output legend**: The output file's *family-name.clinvar-hgmd.binary* and the *family-name.gnomad.binary* have the following structure:
  1. *Index*: Absolute position of the aligned aminoacids.
  2. *Gene-Sequence*: Canonical protein sequence in the format: “aminoacid_position” for each of the genes belonging to the family (one per column).
  4. *Parazscore*: Paralog score observed for that Index position.
  5. *Gene-Binary Annotation*: At least one missense variant observed at corresponding aminoacid (YES=1, NO=0). One column per gene-family member.
  6. *Gene:Disease*: Gene ID coupled with the corresponding disease association observed. This is a collapsed field and more than one gene can have disease associations at the same index positions. Since no disease associations are present in the *family-name.gnomad.binary* file, the tag "GeneID:gnomad" is collapsed on this column. 
### Part 2:
- **Command**: From terminal and inside the repository directory run:

  `$ Rscript Part-2-burden-analysis.R`

- **Output**: The R script will calculate the missense burden analysis over the whole alignment and identify PERs when the difference between the burden of the general population’s missense variants and patient’s pathogenic missense variants becomes significant. Burdens plots are produced in the same format as the one shown in the PER viewer (http://per.broadinstitute.org/). The complete burden analysis are written in a single file with “bin9.stats” extension denoting the 9 amino acid bin size used for the burden analysis. 
- **Output legend**: The file *family-name.bin9.stats* has the following structure:
  1. *Index*: Absolute position of the aligned aminoacids.
  2. *Sequence*: Canonical protein sequence in the format: “aminoacid_position” for each of the genes belonging to the family (one per column).
  3. *Parazscore*: Paralog score observed for that Index position.
  4. *Adj_bin_count*: Adjusted burden of missense variants observed in the general population.
  5. *DM_adj_bin_count*: Adjusted burden of missense variants observed in patients from CLinvar-HGMD (Pathogenic, Likely Pathogenic and/or Disease Mutations).
  7. *Gene:Disease*: Gene ID coupled with the corresponding disease association observed. This is a collapsed field and more than one gene can have disease associations at the same index positions.  
  9. *fisher.p*: Nominal p-value from fisher exact test. 
  10. *or*: odd ratio.
  11. *ci1*: 95% lower confidence interval.
  12. *ci2*: 95% upper confidence interval.
  13. *adj.p*: Bonferroni adjusted p-value considering the amount of bins tested in the alignment (the longer the alignment, lower the alpha).
  14. *log.adj.p*: Logarithm of the adjusted p-value.
  15. *aa.per*: Index position belongs to a PER.
  16. *proxy*: Index position is the anchor of the bin tested. The bin size is 9 aminoacids, with the structure -4aa anchor-aa +4aa. Anchor aminoacid determine the stats of the whole bin. 
  18. *per.tag*: Number of PER. If overlapping bins are significant, the bins are fused together keeping the strongest proxy. 
  19. *per.start*: Start index Position of the PER.
  20. *per.end*: End index position of the PER.
  21. *size*: Aminoacid size of corresponding PER.
  
 ## Final remarks
After running the provided example the user should be able to produce the files contained in "*Example-output*" folder. The complete study and detailed method description is currently available as a preprint entitled: “Identification of pathogenic variant enriched regions across genes and gene families” (https://www.biorxiv.org/content/10.1101/641043v1).
