cellranger="/home/ye/Software/cellranger/cellranger"
output="/home/ye/Data/Zoc/Cell/10X/TUTORIAL/data/fastq/tiny-bcl"
transcriptome="/home/ye/Data/Zoc/Cell/10X/TUTORIAL/refdata-cellranger-GRCh38-3.0.0"

$cellranger count --id test_sample \
	--fastqs $output \
	--sample test_sample \
	--transcriptome $transcriptome \
	--chemistry "threeprime"  # 3'
