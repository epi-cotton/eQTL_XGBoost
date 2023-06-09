

### 1 Identification of Single Nucleotide Polymorphisms (SNPs)

The genome reference is TM-1 v2.1 and can be obtained from the source ( http://cotton.zju.edu.cn/source/TM-1_V2.1.fa.gz)

#####1.1 The whole genome DNA sequences (WGS) were quality trimmed using fastp (version 0.12.2)

```shell
fastp -i WGS_read_1.gz -I WGS_read_2.gz -o WGS_read_qc_1.gz -O WGS_read_qc_2.gz
```

##### 1.2 The quality summary of reads was reported using FastqCount_v0.5.（version 0.5）

```shell
FastqCount_v0.5 WGS_read_1.gz > WGS_read_report.txt
```

##### 1.3 The genome index was built using BWA for DNA alignment (version 0.12.2)

```shell
bwa index ref.fa ref.fa -a bwtsw
```

##### 1.4 DNA Reads were aligned to the reference using bwa(version 0.12.2),  samtools (Version: 1.6 (using htslib 1.6)) and picard (version 1.124) 

```shell
destdir=./bam
for fname in *_1.fastq.gz
do
base=${fname%_1.fastq.gz}
sm=${fname%%_*}
bwa mem -t 18 -M -R "@RG\tID:${sm}\tSM:${sm}\tPL:illumina\tLB:${sm}" ref.fa "$fname" "${base}_2.fastq.gz" | gzip -3 > "$destdir/${base}.sam.gz"
samtools fixmate -O bam $destdir/${base}.sam.gz $destdir/${base}_fixmate.bam
rm $destdir/${base}.sam.gz
samtools sort -@ 8 -O bam -o $destdir/${sm}_sorted.bam -T $destdir/${sm}_temp  $destdir/${base}_fixmate.bam
rm $destdir/${base}_fixmate.bam
java -jar picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT="$destdir/${sm}_sorted.bam" OUTPUT="$destdir/${sm}_dedup.bam" METRICS_FILE="$destdir/${sm}_metrics.txt"
samtools index $destdir/${sm}_dedup.bam
rm $destdir/${sm}_sorted.bam
done
```

##### 1.5 Call SNP using samtools (Version: 1.6 (using htslib 1.6)) and bcftools Version: 1.8 (using htslib 1.8)

The name of 279 samples were sorted in filename.txt 

The BAM files generated in step 1.5 were in the same directory.

```shell
samtools mpileup -q 20 -Q 15 -ugf ref.fa -b filename.txt|bcftools call -vmO z -o raw.vcf.gz
```

#### 1.6 SNP filitering using beagle (version r1399) and vcftools ( version 0.1.16)

```shell
vcftools --gzvcf raw.vcf.gz --remove-indels --recode --recode-INFO-all --out raw_01.vcf
awk '$5 !~ /([[:alpha:]])+,[[:alpha:]]/{print}' raw_01.vcf > raw_02.vcf
vcftools --vcf raw_02.vcf --max-missing 0.8 --minDP 3 --maf 0.05 --recode --recode-INFO-all --out raw_03.vcf
java -jar ~/bin/beagle.r1399.jar gt=raw_03.vcf out=qc.vcf nthreads=5
```

### 2 LncRNA annotation

2.1 The RNA-seq samples were quality-trimmed using fastp (version 0.12.2)

```shell
fastp -i RNA_read_1.gz -I RNA_read_2.gz -o RNA_read_qc_1.gz -O RNA_read_qc_2.gz
```

2.2 The quality summary of reads were reported using FastqCount_v0.5 （version 0.5）

```shell
FastqCount_v0.5 RNA_read_1.gz > RNA_read_report.txt
```

2.3 The genome reference index were build using hisat2（version 2.1.0）

```shell
gffread gene.gff3 -T -o gene.gtf
hisat2_extract_splice_sites.py gene.gtf > gene.ss
hisat2_extract_exons.py gene.gtf > gene.exon
hisat2-build -p 20 --ss gene.ss --exon gene.exon genome.fa genome.fa
```

2.4 LncRNA assemble using cufflinks (verion 2.2.1) 

The clean read of transcriptome were aligned to the reference genome. 

The mapping ratio were saved in the sample_summray.txt.

```shell
hisat2 -p 8 --dta -x ref.fa -1 RNA_read_qc_1.gz -2 RNA_read_qc_2.gz -S sample.sam --summary-file sample_summray.txt
samtools sort -@ 5 -o sample.sorted.bam sample.sam

```

lncRNA assemble 

The novel transcript were assembled using StringTie ( v1.3.3b); 

```shell
ls *bam|while read sample; do stringtie -G gene.gtf -p 10 -o ${sample}.gtf $sample;done
ls *gtf > mergelist.txt 
stringtie --merge -o merged.gtf -c 3 ./mergelist.txt 
cuffcompare -r gene.gtf -p 4 merged.gtf -o merged_lncRNA
awk '$3 == "x"|| $3 == "u"|| $3 == "i" {print $0}'- merged_lncRNA.merged.gtf.tmap > novel.gtf.tmap 
```

Coding potential evaluation using CPC v2.0 ，pfam_scan (version 1.6.2)

LncRNA with length > 200 bp 

```shell
#awk '$3 == "u" {print $0}' merged_lncRNA.merged.gtf.tmap > novel.gtf.tmap 
awk '$10 >200 {print $0}' novel.gtf.tmap > novel.longRNA.gtf.tmap 
awk '{print $5}' novel.longRNA.gtf.tmap | perl extract_gtf_by_name.pl merged.gtf - > novel.longRNA.gtf 
gffread -g genome.fa -w lncRNA.fa ./novel.longRNA.gtf 
TransDecoder.LongOrfs -t lncRNA.fa 
# This step generated a file named longest_orfs.pep
pfam_scan.pl -cpu 8 -fasta ./exon.fa.transdecoder_dir/longest_orfs.pep -dir PfamScan/base/ > pfam_scan_lncRNA.txt
perl -ne 'print if /Domain/' pfam_scan_lncRNA.txt |perl -ne 'print $1."\n" if /::(.*?)::/' > Transcript_with_known_domain.txt
CPC2.py -i exon.fa -o cpc_output.txt
perl -ne 'print if /noncoding/' cpc_output.txt |cut -f 1 > Transcript_with_no_coding_potential.txt
cat Transcript_with_no_coding_potential.txt Transcript_with_known_domain.txt |sort|uniq -d > temp.txt
sort Transcript_with_no_coding_potential.txt temp.txt temp.txt |uniq -u > lncRNA_list.txt
cat lncRNA_list.txt| perl extract_gtf_by_name.pl merged.gtf - > LncRNA.gtf
cat  LncRNA.gtf ref.gtf > lncRNA_mRNA.gtf

```

#### supplementary script 1: extract_gtf_by_name.pl

```shell
#!/usr/bin/perl -w
use strict;
my $gtf=shift;
open FH,$gtf or die;
my %map;
while(<FH>){
	chomp;
	next if /^#/;
	my @field=split "\t";
	#print $field[8]."\n";
	$field[8]=~/transcript_id \"(.+?)\"/;
	my $tid=$1;
	if(!exists $map{$tid}){
		$map{$tid}=$_;
	}else{
		$map{$tid}.="\n".$_;
	}
}
my $idfile=shift;
open FH,$idfile or die;
while(<FH>){
	chomp;
	print $map{$_}."\n";
}
```

##### 2.5  Estimate transcript abundances

```shell
ls *bam|while read id; do stringtie -A ${id}_fpkm.txt -p 8 -G lncRNA_mRNA.gtf -o temp.gtf $id;done ls *fpkm.txt|while read id; do perl -ne 'print unless /^STRG/' $id > 1 && mv 1 $id;done
```

#### 3 eQTL analysis

The gene expression values were formatted in a matrix named gene_expression.matrix.txt.

#### 3.1 Identification of genes expressed in more than 5% accession

```shell
perl ./gene_expression_filter.pl gene_expression.matrix.txt > QC_gene_expression.matrix.txt 
```

supplementary script 2: gene_expression_filter.pl

```shell
open (IN,"$ARGV[0]"); # read fpkm.txt
$head= <IN>; 
print $head;
$sample_num=0;
$sample_num_expressed=0;
while (<IN>) {
chomp;
@a=split("\t",$_);
$id=shift@a;
foreach $element(@a){
	if ($element > 1) { $sample_num_expressed++;}
	    $sample_num ++;
}
	$raio= $sample_num_expressed/$sample_num;
	if ($raio > 0.05) { print $_."\n";
	}

$sample_num=0;
$sample_expressed=0;
```

#### 3.2 Normalization of gene expression values

```
Rscrpit ./qq_normal.r gene_expression.matrix.txt > QC_gene_expression.matrix.txt 
```

supplementary script 3: qq_normal.r

```r
setwd("./")
df=read.table("QC_gene_expression.matrix.txt ",header = T,sep = "\t",check.names = F,row.names = 1)
df_2=df
for (i in 1:nrow(df)) {
  df_2[i,] = qqnorm(df[i,])$x
}
write.table(df_2,"QC_gene_expression.matrix.txt 2_qqnorm.txt",quote = F,sep = "\t")
```

#### 3.3 eQTL using emmax https://genome.sph.umich.edu/wiki/EMMAX 

The gene expression of each gene was normalized, and the import format was changed to EMMAX.

```perl
perl ./emmax_import.pl QC_gene_expression.matrix.txt 
```

```shell
open (IN,"$ARGV[0]");
$head=<IN>;
chomp $head;
@head=split("\t",$head);
while ($line=<IN>) {
chomp $line;
@b=split("\t",$line);
$id=shift@b;
$i=0;
open (OUT,">${id}_time.txt");
foreach $element (@b) {
print OUT $head[$i]."\t";
print OUT $head[$i]."\t";
print OUT $element."\n";
$i++
}
}
```

Association analysis were conducted using EMMAX

```shell
plink --vcf DNA.vcf --recode --out eQTL --noweb
plink --file eQTL --recode12 --output-missing-genotype 0 --transpose --out eQTL --noweb
emmax-kin -v -d 10 eQTL # this step generated a file named eQTL.BN.kinf
plink --noweb --file eQTL --make-bed --out eQTL
gcta64 --bfile eQTL --make-grm --autosome --out eQTL
gcta64 --grm eQTL --pca --out file_temp_pca
$cut -d " " -f 2,3,4,5 file_temp_pca.eigenvec|awk '{print $1"\t"$1"\t1\t"$2"\t"$3"\t"$4}' > cov.txt
vcftools --vcf DNA.vcf --missing-indv ## this step generated a out.imss file
cut -f 1 out.imss |grep -v "INDV" > id_order.txt
# Ensure that the ID is consistent between genotype and phenotypic data.
ls phenotype*txt|while read id; do csvtk join -t -H ../id_order.txt $id > 1 && mv 1 $id;done
# Assocation analysis for each gene
ls *txt |while read id; do emmax -v -d 10 -t eQTL -p $id -k eQTL.BN.kinf -o $id -c ../cov.txt";done
## this step will generated GWAS result for each gene
```

#### 3.4 Calculating the best Pvalue 

```shell
vcftools --vcf DNA.vcf --plink --out output
plink --noweb --file input --dog --transpose --recode12 --out output
plink --noweb --file input --dog --make-bed --out output
java -Xmx4g -jar ~/biosoftware/gec/gec.jar --no-web --effect-number --plink-binary input --genome --out output
# The best pvalue of this study is 2.18E-6
```

#### 3.5 Manhandun Plot for eQTL

The result of emmax were used for Manhandun Plot using R package QQman

```R
#The output was generated from emmax
# Only suit for my cotton project
argv <- commandArgs(TRUE)
df=read.table(argv[1],header = F,sep="\t")
head(df)
library(dplyr)
library(tidyr)
library(qqman)
head(df)
new_df=df %>% separate(V1,c("Chr","Bp"),":",remove=F)
df_gwas =subset(new_df,select=c("V1","Chr","Bp","V3"))
colnames(df_gwas)=c("SNP","CHR","BP","P")
df_gwas$CHR=as.numeric(df_gwas$CHR)
df_gwas$BP=as.numeric(df_gwas$BP)
png(paste(argv[1],"_manhattan.png",sep=""),width= 1000,height = 500)
manhattan(df_gwas,col = c("#e89c39", "#459a70"),suggestiveline = F, genomewideline = F)
dev.off()
png(paste(argv[1],"_qqplot.png",sep=""))
qq(df_gwas$P,col="#e89c39")
dev.off()

```

#### 3.6 Identification of feature SNP 

Identification of feature eSNP 

```shell
mkidr 01raw_data && cd 01raw_data
ls *ps|while read id; do awk '$3 < 2.18E-6 {print $0}' ${id} > ${id}_qc.txt";done
find . -name "*" -type f -size 0c | xargs -n 1 rm -f  
ls *qc.txt|while read id; do perl qc_ps_2_bed.pl $id;done
for i in *bed;do mv $i ../02bed ;done
cd ../02bed
for i in *bed;do mergeBed -i $i -d 100000 -c 1,4 -o count,collapse -delim ";" > ${i}_merge.txt;done
mkdir 03qc_file && cd 03qc_file
for i in *txt; do perl -ne 'print if /;.*?;/' $i > ../03qc_3/$i ;done

for i in *merge.txt; do cut -f 5 $id > 1 && mv 1 $id;done
perl -pi -e 's/;/\n/g' *merge.txt
perl -pi -e 's/,/\t/g' *merge.txt
```

supplementary script 3: qc_ps_2_bed.pl

```perl
open (IN,"$ARGV[0]");
open (OUT,">$ARGV[0].bed");
while (<IN>) {
chomp;
@a=split("\t",$_);
$id=$a[0];
@loc_info=split(":",$id);
$chr=$loc_info[0];
$loc=$loc_info[1];
$pvalue=$a[3];   
print OUT $chr."\t".$loc."\t".$loc."\t".$id.",".$pvalue."\n";
}
```

LD clump

```shell
ls *merge.txt|while read id; do cut -f 1 $id;done |sort|uniq > Uniq_SNPs.txt
plink --noweb --file file --make-bed --out file
plink --noweb --bfile file --extract Uniq_SNPs.txt --make-bed --out eQTL_snp
plink --bfile eQTL_snp --ld-window-kb 1000 --ld-window-r2 0.1 --ld-window 99999 --r2
awk 'BEGIN{OFS="\t"}{$1=$1;print}' plink.ld > 1 && mv 1 plink.ld
cut -f 3,6,7 plink.ld > snp-snp-r2.txt
perl eQTL_LD_clump.pl 
awk '{print FILENAME"\t"$0}' *rm_dup.txt|perl -pi -e 's/.txt.ps_qc.txt.bed_merge.txt_rm_dup.txt//' > eQTL.txt
```

supplementary script 4: eQTL_LD_clump.pl 

```shell
open (IN,"snp-snp-r2.txt");     #读取snp对
while (<IN>) {
chomp;
@a=split("\t",$_);
@arrary_temp= sort($a[0],$a[1]);
$pair= join("-",@arrary_temp);
#print $pair."\n";
$hash_pair{$pair}=1; 
}
## 处理目录下所有文件
my $DIR_PATH="./";
opendir DIR, ${DIR_PATH} or die "Can not open";
@filelist = readdir DIR;
foreach $file (@filelist) {
if ($file !~ /Qc/) {next;}  
open (IN2,"$file"); 
open (OUT,">${file}_rm_dup.txt");
$first= <IN2>;
chomp $first;
@b=split("\t",$first);
$snp_temp= $b[0];
$snp_temp_pvalue= $b[1];
while ($line= <IN2>) {
chomp $line;
@c=split("\t",$line);
$snp= $c[0];
$snp_pvalue= $c[1];
@arrary_temp_2= sort($snp,$snp_temp);
$str= join("-",@arrary_temp_2);
if (exists $hash_pair{$str}  ) { 
if ($snp_pvalue  < $snp_temp_pvalue){
$snp_temp = $snp;
$snp_temp_pvalue= $snp_pvalue;}
}
else { 
print OUT $snp_temp."\t".$snp_temp_pvalue."\n";
$snp_temp = $snp;
$snp_temp_pvalue= $snp_pvalue;
}
}
print OUT $snp_temp."\t".$snp_temp_pvalue."\n";
}#
```

##### 3.7 LD block show of eQTL

```shell
LDBlockShow -InVCF DNA.vcf.gz -OutPut output.png -Region D08:2212434-2429470 -InGWAS ./eQTL.ps -InGFF gene.gff
```

#### 4 GRN construction

The script similar to eQTL_LD_clump.pl 

#### 5 Heritability estimation

Three data set of SNP were generated. The number of SNPs were consists.

```shell
for i in $(seq 1 100)  
do
gcta64 --bfile eQTL --extract snp_dataset1.txt --make-grm --out 01grmResult
gcta64 --bfile eQTL --extract snp_dataset2.txt --make-grm --out 02grmResult
gcta64 --bfile eQTL --extract snp_dataset3.txt --make-grm --out 03grmResult
gcta64 --reml --grm ../01grmResult --pheno ./phenotype/SI07an.txt --out ${i}_SI07an --reml-maxit 
10000 --reml-alg 2
done

```

#### 6 XGBoost for Phenotypic prediction

```python
import os
import math
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import xgboost as xgb
from xgboost import plot_importance
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import collections, functools, operator
  data = pd.read_csv("Input.txt", sep='\t')
trait=  data.loc[:, 'trait']
data=data.iloc[:,2:]
X, y = data, trait
X= np.log(X+0.1)
params = {
    'booster': 'gbtree',
    'objective': 'reg:gamma',
    'gamma': 0.1,
    'max_depth': 5,
    'lambda': 3,
    'subsample': 0.7,
    'colsample_bytree': 0.7,
    'min_child_weight': 3,
    'silent': 1,
    'eta': 0.1,
    'seed': 1000,
    'nthread': 4,
}
params2 = {
        'booster': 'gbtree',
        'objective': 'reg:linear',
        'gamma': 50,
        'max_depth': 6,
        'lambda': 3,
        'subsample': 0.25,
        'colsample_bytree': 0.33,
        'min_child_weight': 1,
        'silent': 1,
        'eta': 0.075,
        'seed': 1000,
        'nthread': 4,
    }

SumScore=[]
for i in range(0,100):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=i+1000)
    seed=1000+i
    
    dtrain = xgb.DMatrix(X_train, y_train)
    num_rounds = 300

    model = xgb.train(params, dtrain,num_rounds )
    dtest = xgb.DMatrix(X_test)
    ans = model.predict(dtest)
    mse = mean_squared_error(ans, y_test)
    r2 = r2_score(y_test,ans,sample_weight=None)

    print('rmse = %.4f' % mse)
    print ('R2 = %.4f' % r2)

    with open("Output_R2.txt", 'a') as f2:
        f2.write('R2 = %.4f' % r2+"\n")
    SumScore.append(model.get_score(importance_type='weight'))

# sum the values with same keys
result = dict(functools.reduce(operator.add,
         map(collections.Counter, SumScore)))
 
with open("Output.txt", 'w') as f: 
    for key, value in result.items(): 
        f.write('%s\t%s\n' % (key, value))

ax = plot_importance(result,max_num_features=20)
plt.tight_layout()
ax.figure.savefig('Out.pdf')

```





