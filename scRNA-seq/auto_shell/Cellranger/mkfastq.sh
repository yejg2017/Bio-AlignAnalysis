run="/home/ye/Data/Zoc/Cell/10X/TUTORIAL/tiny-bcl"
csv="/home/ye/Data/Zoc/Cell/10X/TUTORIAL/cellranger-tiny-bcl-simple-1.2.0.csv"
output="/home/ye/Data/Zoc/Cell/10X/TUTORIAL/data/fastq/tiny-bcl"

bclfastq="/home/ye/Software/Python/location/bin/bcl2fastq"
cellranger="/home/ye/Software/cellranger/cellranger"


$cellranger mkfastq --run $run \
	--csv $csv \
	--id "tiny-bcl" \
	--output-dir=$output
