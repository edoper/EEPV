# EEPV
A framework to identify and visualize Exons Enriched with Pathogenic Variants

## Introduction

Using publicly available patient and general population variants we developed a framework to identify exons significantly enriched with pathogenic variation on a genome-wide level. We detected 179 and 79 genes with an exon enriched for pathogenic and likely pathogenic SNV and INDELs, respectively. All these exons harbor >50% of all pathogenic variation observed within their corresponding gene. Our framework will facilitate efficient candidate gene selection for the development of exon-skipping therapies.

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

Since the GnomAD exome data is divided by chromosomes, we merge them and create a merged vcf file.
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
We repeat the same process with SNVs.
```
# Extracting SNVs.
for S in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
	echo "working with $S"
	bcftools view --types snps gnomad.exomes.r2.1.1.sites.$S.vcf.bgz >$S.snps.exome.vcf
	bgzip -c $S.snps.exome.vcf >$S.snps.exome.vcf.gz
	tabix -p vcf $S.snps.exome.vcf.gz
done
bcftools view --types snps gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz >all.snps.genome.vcf
bgzip -c all.snps.genome.vcf >all.snps.genome.vcf.gz
tabix -p vcf all.snps.genome.vcf.gz

# Concatenating exome SNVs.
bcftools concat -O z -o all.snps.exome.vcf.gz 1.snps.exome.vcf.gz 2.snps.exome.vcf.gz 3.snps.exome.vcf.gz 4.snps.exome.vcf.gz 5.snps.exome.vcf.gz 6.snps.exome.vcf.gz 7.snps.exome.vcf.gz 8.snps.exome.vcf.gz 9.snps.exome.vcf.gz 10.snps.exome.vcf.gz 11.snps.exome.vcf.gz 12.snps.exome.vcf.gz 13.snps.exome.vcf.gz 14.snps.exome.vcf.gz 15.snps.exome.vcf.gz 16.snps.exome.vcf.gz 17.snps.exome.vcf.gz 18.snps.exome.vcf.gz 19.snps.exome.vcf.gz 20.snps.exome.vcf.gz 21.snps.exome.vcf.gz 22.snps.exome.vcf.gz X.snps.exome.vcf.gz Y.snps.exome.vcf.gz
tabix -p vcf all.snps.exome.vcf.gz

# PASS filtering.
bcftools view -O z -f 'PASS' -o all.snps.exome.vcf.pass.gz all.snps.exome.vcf.gz
bcftools view -O z -f 'PASS' -o all.snps.genome.vcf.pass.gz all.snps.genome.vcf.gz
tabix -p vcf all.snps.exome.vcf.pass.gz
tabix -p vcf all.snps.genome.vcf.pass.gz

# Vcf merging.
bcftools merge -O z -o snps.vcf.gz all.snps.exome.vcf.pass.gz all.snps.genome.vcf.pass.gz
```
After obtaining the indels and SNVs vcf files, we will use an auxiliary perl script `cleaning-gnomaAD.pl` to extract specific information from the vcf files. The final output is a bed file that will be used later.
```
# Gunzip vcfs.
# Extracting information from indels vcf.
perl cleaning-gnomaAD.pl indels.vcf indels.gnomad.bed
# Reformat snps, split by chromosome, and extract information.
zcat snps.vcf.gz | awk '!/^#/{print>$1}'
for F in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
	echo "working with $F"
	perl cleaning-gnomaAD.pl $F $F.bed
done
# Merge snps.
cat 1.bed 2.bed 3.bed 4.bed 5.bed 6.bed 7.bed 8.bed 9.bed 10.bed 11.bed 12.bed 13.bed 14.bed 15.bed 16.bed 17.bed 18.bed 19.bed 20.bed 21.bed 22.bed X.bed Y.bed >cat.snps.bed
grep -v -P '^chr' cat.snps.bed >snv.gnomad.bed
# Format change.
awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6":"$7}' snv.gnomad.bed >snv.gnomad.bed
awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5":"$6":"$7}' indels.gnomad.bed >indel.gnomad.bed   
```
### Part 2: Pathogenic variants from ClinVar and HGMD.
The second part requires datasets downloaded from ClinVar `variant_summary.txt` (https://www.ncbi.nlm.nih.gov/clinvar/) and HGMD `hg19_hgmd.txt` (http://www.hgmd.cf.ac.uk/ac/index.php). We use the auxiliary scripts `source2input-hgmd.pl` and `source2input-clinvar.pl` to generate bed files with SNVs and indels.
```
# Extracting data
perl source2input-hgmd.pl
perl source2input-clinvar.pl 
perl -pi -e 's/\;/\|/g' snv.clinvar.bed
perl -pi -e 's/\;/\|/g' indel.clinvar.bed
perl -pi -e 's/\?/o/g' snv.clinvar.bed
perl -pi -e 's/\?/o/g' indel.clinvar.bed
# Change chromosomes into a numeric format
perl -pi -e 's/^X\t/23\t/g' snv.clinvar.bed
perl -pi -e 's/^X\t/23\t/g' indel.clinvar.bed
perl -pi -e 's/^X\t/23\t/g' snv.hgmd.bed
perl -pi -e 's/^X\t/23\t/g' indel.hgmd.bed
perl -pi -e 's/^Y\t/24\t/g' snv.clinvar.bed
perl -pi -e 's/^Y\t/24\t/g' indel.clinvar.bed
perl -pi -e 's/^Y\t/24\t/g' snv.hgmd.bed
perl -pi -e 's/^Y\t/24\t/g' indel.hgmd.bed
# Merge ClinVar and HGMD datasets.
cat snv.clinvar.bed snv.hgmd.bed >snv.patogenic.bed
cat indel.clinvar.bed indel.hgmd.bed >indel.patogenic.bed
```
Here we use an auxiliar script called `multiple-column-sort.pl`, the output has the same name as the input with a `.sorted` suffix.
```
perl multiple-column-sort.pl snv.patogenic.bed 2n 1n
perl multiple-column-sort.pl indel.patogenic.bed 2n 1n
# Removing MT and reformat.
grep -v -P '(^MT)|(\tna\tna\t)' snv.patogenic.bed.sorted >snv.patogenic.ungrouped
grep -v -P '(^MT)|(\tna\tna\t)' indel.patogenic.bed.sorted >indel.patogenic.ungrouped
perl -pi -e 's/\,/_THIS_IS_A_COMMA_/g' snv.patogenic.ungrouped
perl -pi -e 's/\,/_THIS_IS_A_COMMA_/g' indel.patogenic.ungrouped
# Group to get unique rows.
bedtools groupby -i snv.patogenic.ungrouped -g 1,2,3 -c 4,5 -o collapse,collapse >input.patogenic.snv
bedtools groupby -i indel.patogenic.ungrouped -g 1,2,3 -c 4,5 -o collapse,collapse >input.patogenic.indel
bedtools groupby -i snv.gnomad.bed -g 1,2,3 -c 4,5 -o collapse,collapse >input.gnomad.snv
bedtools groupby -i indel.gnomad.bed -g 1,2,3 -c 4,5 -o collapse,collapse >input.gnomad.indel
# Replace commas.
perl -pi -e 's/\,/;/g' input.patogenic.snv
perl -pi -e 's/\,/;/g' input.patogenic.indel
perl -pi -e 's/\,/;/g' input.gnomad.snv
perl -pi -e 's/\,/;/g' input.gnomad.indel
# Reformat with original commas and chromosomes from phenotypes.
perl -pi -e 's/_THIS_IS_A_COMMA_/\,/g' input.patogenic.snv
perl -pi -e 's/_THIS_IS_A_COMMA_/\,/g' input.patogenic.indel
#perl -pi -e 's/\A23\t/\AX\t/g' input.patogenic.snv
#perl -pi -e 's/\A24\t/\AY\t/g' input.patogenic.indel
```
### Part 3: Transcript mapping pathogenic + gnomAD
Exon data was downloaded from UCSC. We provide two examples, the APP gene located on chr21 and the GABRB2 gene located on chr5, on a custom `*.map` format.
```
for C in 5 21; do
	echo "working with $C"
	awk -F"\t" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' chr$C.map >$C.backbone
done
```
From the `*.map` file we generate a `*.backbone` file, and we use two auxiliary scripts, one for SNVs, and another for indels.
```
perl mapping.snvs.pl
mkdir map-snv
for C in 5 21; do
	mv ${C}.map ./map-snv/
done

perl mapping.indel.pl
mkdir map-indel
for C in 5 21; do
	mv ${C}.map ./map-indel/
done
```
We have created two different folders, one for SNVs and another for indels. Now we generate `*.map` files for each chromosome in the two created folders, group the variants in each, and count the variants in each exon.
```
for d in map-snv map-indel; do
cd ${d}
for I in 5 21; do 
	echo "chr $I"
	perl -pi -e 's/,/#/g' $I.map
	#bedtools groupby -i chr$I.map -g 1,2,3,10,9,8 -c 4,4,4,16,21,22,24 -o count,min,max,sum,sum,sum,collapse >chr$I.counts
	bedtools groupby -i $I.map -g 1,2,3,10,9,8 -c 4,4,4,16,16,20,22 -o count,min,max,sum,sum,sum,collapse >chr$I.counts
	perl -pi -e 's/\|/;/g' chr$I.counts;
	perl -pi -e 's/,/;/g' chr$I.counts;
	perl -pi -e 's/#/,/g' chr$I.counts;
	perl -pi -e 's/0;//g' chr$I.counts;
	perl -pi -e 's/;0\n/\n/g' chr$I.counts;
	perl -pi -e 's/\t;/\t/g' chr$I.counts;
	awk -F"\t" '{print $1"\t"$2}' chr$I.counts | sort | uniq -c | awk '{print $1"\t"$2}' >chr$I.genes
	perl -pi -e 's/#/,/g' $I.map;
	echo "done with chr$I";
done
cd ..
done
```

### Part 4: Exome enrichment testing
Here we check for variant enrichment on each exon. We check for SNVs enrichment first, and indel enrichment second. We use two auxiliary scripts, `counting.by.exons.pl` and `testing.exons.R`.
```
# SNVs.
cd map-snv
cat *.counts >all-exon.txt
cat *.genes | awk -F"\t" '{print $2}' >all-genes.txt
cp ../counting.by.exons.pl ./
cp ../testing.exons.R ./
perl counting.by.exons.pl
Rscript testing.exons.R
cd ..
# Indels.
cd map-indel
cat *.counts >all-exon.txt
cat *.genes | awk -F"\t" '{print $2}' >all-genes.txt
cp ../counting.by.exons.pl ./
cp ../testing.exons.R ./
perl counting.by.exons.pl
Rscript testing.exons.R
cd ..
```
We extract the results using the auxiliary script `results.R`
```
Rscript results.R
```

### Part 5: Burden analysis
We want to analyze the effect of each variant on each exon. Here we use the auxiliary script `burden-analysis1.R` which does the analysis.
```
for i in map-snv map-indel; do
cd ${i}
for ii in 5 21; do
mkdir chr${ii}
cp ${ii}.map ./chr${ii}
done
cp ../burden-analysis1.R ./
Rscript burden-analysis1.R
cd ..
done
```
### Part 6: Plotting results
The final part is plotting the results.
```
for i in map-snv map-indel; do
cd ${i}
mkdir all-gene-burden-main-results
mkdir all-genes-burden-introns-input
mv ./chr5/*.burden ./all-gene-burden-main-results/
mv ./chr21/*.burden ./all-gene-burden-main-results/
cp ../add.introns.pl ./all-gene-burden-main-results/
cd all-gene-burden-main-results
perl add.introns.pl
cd ..
cd ..
done

# Esto escribe resultados en all-genes-burden-introns-input
# In snv
cp toplot-snv.pl ./map-snv/all-genes-burden-introns-input
cd ./map-snv/all-genes-burden-introns-input
mkdir plot-snv
perl toplot-snv.pl
cp *.input ./plot-snv
cp -r plot-snv ../..
cd ../..
# In indel
cp toplot-indel.pl ./map-indel/all-genes-burden-introns-input
cd ./map-indel/all-genes-burden-introns-input
mkdir plot-indel
perl toplot-indel.pl
cp *.input ./plot-indel
cp -r plot-indel ../..
cd ../..
# Esto escribe resultados en plot-snv y plot-indel

for i in snv indel; do
cd plot-${i}
mv 5.input GABRB2.input
mv 21.input APP.input
cd ..
done
```
After processing the data, we use the final r script `plots.R` to generate PDF plots for APP and GABRB2.
```
Rscript plots.R
```

## Final remarks
Interpretation of the results here.
