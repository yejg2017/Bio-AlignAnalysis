# Bio-Mapping
Use the major mapping tools to analysize mouse dataset,like STAR,Bowtie2,samtools...


### 主要是使用 STAR，Bowtie等这些工具来做mapping
* 1) [auto_py](./auto_py)目录下的scripts主要是使用python来生成shell代码，实现并行的mapping(然而主要问题是很消耗内存）
* 2) [auto_shell](./auto_shell)目录下的scripts主要是实现STAR，Bowties2，Cellranger等一些常用工具来mapping实现的demo过程

### Attention
[reference](https://asia.ensembl.org/info/data/ftp/index.html)主要的*gtf*和*fasta*等文件都需要提前在这网站先下载
* 1) ftp://ftp.ensembl.org/pub/release-95/
* 2) https://asia.ensembl.org/info/data/ftp/index.htm


这些软件的具体是使用，自己参考吧，google去

* [auto_shell](./auto_shell) 里面的directory都是个人实现create的，比如，/home/ye/Data/Zoc/Cell/reference/data/qc/star/align  主要是
[align](auto_shell/STAR/alignReads.sh)之后得到的比对sam文件，主要包括的子目录有B1...,B2...,B3..../Align.out.sam等
