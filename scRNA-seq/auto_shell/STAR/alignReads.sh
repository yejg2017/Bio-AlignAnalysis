root=$1  #  fastq root
star=$2
genomeDir=$3  # /home/ye/Data/Zoc/Cell/reference/data/qc/star/genomeDir 
gtf=$4  #  /home/ye/Data/Zoc/Cell/reference/dna/mouse/gtf/Mus_musculus.GRCm38.95.gtf
output=$5    # output file : /home/ye/Data/Zoc/Cell/reference/data/qc/star/align

if [ -f $output ];
then 
	echo "$output exists"
else
	mkdir -p $output
fi


for i in ` ls $root | grep -E  "*.R1.*fastq.gz"`;
do 
	#echo  $i, ${i%.R1.clean.fastq.gz}_R2.clean.fastq.gz
	if [ ! -f $output/${i%.R1.clean.fastq.gz} ];
	then
		mkdir -p $output/${i%.R1.clean.fastq.gz}
	else
		echo ""
	fi

	$star --runMode alignReads \
		--genomeDir $genomeDir \
		--sjdbGTFfile $gtf \
		--readFilesType Fastx \
		--runThreadN 2 \
		--outTmpKeep None \
		--outSAMtype SAM \
		--outFileNamePrefix  $output/${i%.R1.clean.fastq.gz}/  \
		--outFilterMismatchNmax 100 \
		--seedSearchStartLmax 20 \
		--seedPerReadNmax 100000 \
		--seedPerWindowNmax 1000 \
		--alignTranscriptsPerReadNmax 100000 \
		--alignTranscriptsPerWindowNmax 10000 \
		--readFilesCommand zcat \
		--genomeLoad NoSharedMemory  \
		--readFilesIn $root/$i $root/${i%.R1.clean.fastq.gz}.R2.clean.fastq.gz

done

