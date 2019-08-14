* 1. 首先安装软件

安装conda，使用清华镜像的conda，说明书地址：https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/

下载并安装miniconda，地址：https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/ 

* 2. 更改镜像配置
```shell
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh 
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --set show_channel_urls yes
```

* 3. 然后就可以根据流程来使用conda安装一系列软件
```shell
conda  create -n wes  python=2 bwa
conda info --envs
source activate wes
```
* 4. 下载安装软件之前先搜索该软件名是否正确，地址：https://bioconda.github.io/recipes.html可以用search先进行检索
```shell
conda search sratools
# 保证所有的软件都是安装在 wes 这个环境下面
conda install sra-tools
conda install samtools
conda install -y bcftools vcftools  snpeff
conda install -y multiqc qualimap

# https://software.broadinstitute.org/gatk/download/
https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
cd ~/biosoft
mkdir -p gatk4 &&  cd gatk4
wget  https://github.com/broadinstitute/gatk/releases/download/4.0.6.0/gatk-4.0.6.0.zip
unzip gatk-4.0.6.0.zip
```
* 5. 然后下载100G的必备数据才能使用GATK,熟悉参考基因组及必备数据库
```shell
/public/biosoft/GATK/resources/bundle/hg38/bwa_index/
|-- [ 20K]  gatk_hg38.amb
|-- [445K]  gatk_hg38.ann
|-- [3.0G]  gatk_hg38.bwt
|-- [767M]  gatk_hg38.pac
|-- [1.5G]  gatk_hg38.sa
|-- [6.2K]  hg38.bwa_index.log 
|-- [ 566]  run.sh

/public/biosoft/GATK/resources/bundle/hg38/
|-- [1.8G]  1000G_phase1.snps.high_confidence.hg38.vcf.gz
|-- [2.0M]  1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
|-- [3.2G]  dbsnp_146.hg38.vcf.gz
|-- [3.0M]  dbsnp_146.hg38.vcf.gz.tbi
|-- [ 59M]  hapmap_3.3.hg38.vcf.gz
|-- [1.5M]  hapmap_3.3.hg38.vcf.gz.tbi
|-- [568K]  Homo_sapiens_assembly38.dict
|-- [3.0G]  Homo_sapiens_assembly38.fasta
|-- [157K]  Homo_sapiens_assembly38.fasta.fai
|-- [ 20M]  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
|-- [1.4M]  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```
* 6. 第一步是QC
包括使用fasqc和multiqc两个软件查看测序质量，以及使用trim_galore软件进行过滤低质量reads和去除接头。
```shell
mkdir ~/project/boy
wkd=/home/jmzeng/project/boy
mkdir {raw,clean,qc,align,mutation}
cd qc 
find /public/project/clinical/beijing_boy  -name *gz |grep -v '\._'|xargs fastqc -t 10 -o ./
#假设质量很差，就过滤：

### step3: filter the bad quality reads and remove adaptors. 
mkdir $wkd/clean 
cd $wkd/clean

find /public/project/clinical/beijing_boy  -name *gz |grep -v '\._'|grep 1.fastq.gz > 1
find /public/project/clinical/beijing_boy  -name *gz |grep -v '\._'|grep 2.fastq.gz > 2
paste 1 2  > config
### 打开文件 qc.sh ，并且写入内容如下： 
source activate wes

bin_trim_galore=trim_galore
dir=$wkd/clean
cat config  |while read id
do
        arr=(${id})
        fq1=${arr[0]}
        fq2=${arr[1]} 
        echo  $dir  $fq1 $fq2 
nohup $bin_trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir  $fq1 $fq2 & 
done 

source deactivate
读质量较好的测序数据进行比对
先走测试数据
```
## 先提取小的fq
source activate wes
find /public/project/clinical/beijing_boy  -name *gz |grep -v '\._' > fq.txt
cat fq.txt |while read id ;do (zcat $id|head -10000 > $(basename $id ".gz"));done
## 然后一个个小fq文件比对
sample='7E5239'
bwa mem -t 5 -R "@RG\tID:$sample^CSM:$sample\tLB:WGS\tPL:Illumina"  /public/biosoft/GATK/resources/bundle/hg38/bwa_index/gatk_hg38 7E5239.L1_1.fastq 7E5239.L1_2.fastq  | samtools sort -@ 5 -o 7E5239.bam -
sample='7E5240'
bwa mem -t 5 -R "@RG\tID:$sample^CSM:$sample\tLB:WGS\tPL:Illumina"  /public/biosoft/GATK/resources/bundle/hg38/bwa_index/gatk_hg38  7E5240_L1_A001.L1_1.fastq 7E5240_L1_A001.L1_2.fastq  | samtools sort -@ 5 -o 7E5240.bam -
sample='7E5241'
bwa mem -t 5 -R "@RG\tID:$sample^CSM:$sample\tLB:WGS\tPL:Illumina"  /public/biosoft/GATK/resources/bundle/hg38/bwa_index/gatk_hg38 7E5241.L1_1.fastq 7E5241.L1_2.fastq   | samtools sort -@ 5 -o 7E5241.bam -
## 或者循环比对
# 7E5239    7E5239.L1_1.fastq   7E5239.L1_2.fastq
# 7E5240    7E5240_L1_A001.L1_1.fastq   7E5240_L1_A001.L1_2.fastq
# 7E5241    7E5241.L1_1.fastq   7E5241.L1_2.fastq
INDEX=/public/biosoft/GATK/resources/bundle/hg38/bwa_index/gatk_hg38
cat config |while read id
do

arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2  | samtools sort -@ 5 -o $sample.bam - 

done 
然后走正常的数据

ls /home/jmzeng/project/boy/clean/*1.fq.gz >1
ls /home/jmzeng/project/boy/clean/*2.fq.gz >2
cut -d"/" -f 7 1 |cut -d"_" -f 1  > 0
paste 0 1 2 > config 
source activate wes
INDEX=/public/biosoft/GATK/resources/bundle/hg38/bwa_index/gatk_hg38
cat config |while read id
do

arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
echo $sample $fq1 $fq2 
 bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2  | samtools sort -@ 5 -o $sample.bam -   

done 
最简单的找变异流程
如果要理解参数的意思，参考我的直播基因组 【直播】我的基因组25：用bcftools来call variation

ref=/public/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
source activate wes
time samtools mpileup -ugf $ref  *.bam | bcftools call -vmO z -o out.vcf.gz
ls *.bam |xargs -i samtools index {}
# 没有去除PCR重复
去除PCR重复
下面的代码不能运行，但是你需要理解，理解参数和教程为什么会过时

# samtools markdup -r 7E5241.bam 7E5241.rm.bam
# samtools markdup -S 7E5241.bam 7E5241.mk.bam
完善的GATK流程
source activate wes
GATK=/home/jmzeng/biosoft/gatk4/gatk-4.0.6.0/gatk
ref=/public/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
snp=/public/biosoft/GATK/resources/bundle/hg38/dbsnp_146.hg38.vcf.gz
indel=/public/biosoft/GATK/resources/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  

for sample in {7E5239.L1,7E5240,7E5241.L1}
do 
echo $sample  
# Elapsed time: 7.91 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
    -I $sample.bam \
    -O ${sample}_marked.bam \
    -M $sample.metrics \
    1>${sample}_log.mark 2>&1 


## Elapsed time: 13.61 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
    -I ${sample}_marked.bam \
    -O ${sample}_marked_fixed.bam \
    -SO coordinate \
    1>${sample}_log.fix 2>&1 

samtools index ${sample}_marked_fixed.bam

##  17.2 minutes
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
    -R $ref  \
    -I ${sample}_marked_fixed.bam  \
    --known-sites $snp \
    --known-sites $indel \
    -O ${sample}_recal.table \
    1>${sample}_log.recal 2>&1 

$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"   ApplyBQSR \
    -R $ref  \
    -I ${sample}_marked_fixed.bam  \
    -bqsr ${sample}_recal.table \
    -O ${sample}_bqsr.bam \
    1>${sample}_log.ApplyBQSR  2>&1 

## 使用GATK的HaplotypeCaller命令
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
     -R $ref  \
     -I ${sample}_bqsr.bam \
      --dbsnp $snp \
      -O ${sample}_raw.vcf \
      1>${sample}_log.HC 2>&1  

done 
检查感兴趣基因区域内比对和找变异情况
通过IGV可视化来加深自己对这个流程的把握和理解。

chr17   HAVANA  gene    43044295        43170245
3.5G Jul 21 18:01 7E5240.bam
7.1G Jul 21 21:40 7E5240_bqsr.bam
4.7G Jul 21 20:28 7E5240_marked.bam
4.8G Jul 21 20:44 7E5240_marked_fixed.bam
把这些bam里面的BRCA1基因的reads拿出来：

samtools  view -h  7E5240.bam chr17:43044295-43170245 |samtools sort -o  7E5240.brca1.bam -
samtools  view -h  7E5240_bqsr.bam chr17:43044295-43170245 |samtools sort -o  7E5240_bqsr.brca1.bam -
samtools  view -h  7E5240_marked.bam chr17:43044295-43170245 |samtools sort -o  7E5240_marked.brca1.bam -
samtools  view -h  7E5240_marked_fixed.bam chr17:43044295-43170245 |samtools sort -o  7E5240_marked_fixed.brca1.bam -
ls  *brca1.bam|xargs -i samtools index {}
有了这些特定基因区域的bam，就可以针对特定基因找变异

source activate wes
ref=/public/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
samtools mpileup -ugf $ref   7E5240_bqsr.brca1.bam   | bcftools call -vmO z -o 7E5240_bqsr.vcf.gz
所有样本走samtools mpileup 和bcftools call  流程
仍然是参考我的直播基因组 【直播】我的基因组25：用bcftools来call variation

source activate wes
wkd=/home/jmzeng/project/boy
cd $wkd/align 
ls  *_bqsr.bam  |xargs -i samtools index {}
ref=/public/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
nohup samtools mpileup -ugf $ref   *_bqsr.bam | bcftools call -vmO z -o all_bqsr.vcf.gz & 
比对及找变异结果的质控
raw和clean的fastq文件都需要使用fastqc质控。

比对的各个阶段的bam文件都可以质控。

使用qualimap对wes的比对bam文件总结测序深度及覆盖度。

source activate wes
wkd=/home/jmzeng/project/boy
cd $wkd/clean
ls *.gz |xargs fastqc -t 10 -o ./

cd $wkd/align
rm *_marked*.bam
ls  *.bam  |xargs -i samtools index {} 
ls  *.bam  | while read id ;do (samtools flagstat $id > $(basename $id ".bam").stat);done

conda install bedtools
cat /public/annotation/CCDS/human/CCDS.20160908.txt  |perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}'|sort -u |bedtools sort -i |awk '{print "chr"$0"\t0\t+"}'  > $wkd/align/hg38.exon.bed 

exon_bed=hg38.exon.bed 
ls  *_bqsr.bam | while read id;
do
sample=${id%%.*}
echo $sample
qualimap bamqc --java-mem-size=20G -gff $exon_bed -bam $id  & 
done 

### featureCounts 
gtf=/public/reference/gtf/gencode/gencode.v25.annotation.gtf.gz
featureCounts -T 5 -p -t exon -g gene_id    \
-a  $gtf   *_bqsr.bam -o  all.id.txt  1>counts.id.log 2>&1 & 
比较两个找变异工具的区别
chr1    139213  .       A       G       388     .   
DP=984;VDB=0.831898;SGB=-223.781;RPB=0.599582;MQB=0.0001971;MQSB=0.00457372;BQB=2.59396e-09;MQ0F=0.281504;ICB=0.333333;HOB=0.5;AC=3;AN=6;DP4=371,169,160,84;MQ=21
GT:PL   0/1:198,0,255   0/1:176,0,255   0/1:50,0,158
chr1    139213  rs370723703 A   G   3945.77 .   AC=1;AF=0.500;AN=2;BaseQRankSum=-2.999;ClippingRankSum=0.000;DB;DP=244;ExcessHet=3.0103;FS=2.256;MLEAC=1;MLEAF=0.500;MQ=29.33;MQRankSum=-0.929;QD=16.17;ReadPosRankSum=1.462;SOR=0.863  GT:AD:DP:GQ:PL  0/1:136,108:244:99:3974,0,6459
chr1    139213  rs370723703 A   G   2261.77 .   AC=1;AF=0.500;AN=2;BaseQRankSum=-1.191;ClippingRankSum=0.000;DB;DP=192;ExcessHet=3.0103;FS=9.094;MLEAC=1;MLEAF=0.500;MQ=32.03;MQRankSum=-0.533;QD=11.78;ReadPosRankSum=0.916;SOR=0.321  GT:AD:DP:GQ:PL  0/1:126,66:192:99:2290,0,7128
chr1    139213  rs370723703 A   G   2445.77 .   AC=1;AF=0.500;AN=2;BaseQRankSum=-2.495;ClippingRankSum=0.000;DB;DP=223;ExcessHet=3.0103;FS=10.346;MLEAC=1;MLEAF=0.500;MQ=30.18;MQRankSum=0.486;QD=10.97;ReadPosRankSum=-0.808;SOR=0.300 GT:AD:DP:GQ:PL  0/1:152,71:223:99:2474,0,7901
补充作业
使用freebayes和varscan找变异。

VCF下游分析
主要是：注释和过滤

注释
常用的软件是：VEP，snpEFF，ANNOVAR，我在生信菜鸟团都写过，大家可以根据下面的标题自行检索，去浏览学习

1.Annovar使用记录

2.用annovar对snp进行注释

3.对感兴趣的基因call variation

4.WES（六）用annovar注释

cd ~/biosoft/
ln -s /public/biosoft/ANNOVAR/ ANNOVAR
source activate wes
wkd=/home/jmzeng/project/boy
cd $wkd/mutation 
 mv ../align/*.vcf ./

 ~/biosoft/ANNOVAR/annovar/convert2annovar.pl  -format vcf4old    7E5240_raw.vcf  >7E5240.annovar

~/biosoft/ANNOVAR/annovar/annotate_variation.pl \
-buildver hg38  \
--outfile 7E5240.anno \
7E5240.annovar \
~/biosoft/ANNOVAR/annovar/humandb/
其实并不一定需要vcf文件，软件需要的只是染色体加上坐标即可，对于我们的vcf格式的变异文件， 软件通常会进行一定程度的格式化之后再进行注释。这里的注释主要有三种方式，分别是：

基于基因的注释,exonic,splicing,ncRNA,UTR5,UTR3,intronic,upstream,downstream,intergenic

基于区域的注释, cytoBand,TFBS,SV,bed,GWAS,ENCODE,enhancers, repressors, promoters

基于数据库的过滤,dbSNP,ExAC,ESP6500,cosmic,gnomad,1000genomes,clinvar

补充作业
使用 Variant Effect Predictor 对所有遗传变异进行注释。过滤掉 dbSNP 数据库和千人基因组计划数据库中已知的 SNP。

应用 OMIM 数据库(http://omim.org/)查询蛋白的结构及功能。利用 SIFT ，PolyPhen-2 以及 PROVEAN 软件, 预测 SNV 对蛋白质功能的影响程度，仅当 3 种软件均预测同一遗传变异对蛋白质的功能影响较大时，才认定该遗传变异具有高危害性。利用 PROVEAN 软件 预测 Indel 对蛋白质功能的影响。

其实dbNSFP数据库，就注释了这些变异对蛋白功能的影响。

变异位点的过滤
使用 GATK的 Joint Calling过滤

通过搜索我们找到了gatk3代码：http://baserecalibrator1.rssing.com/chan-10751514/all_p13.html，但事实上这样的教程过时了，可以抛弃
