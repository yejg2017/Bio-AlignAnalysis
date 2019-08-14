/home/ye/Software/Python/location/envs/Align/bin/STAR --runMode alignReads \
	--genomeDir /home/ye/Data/Zoc/Cell/reference/data/qc/star/genomeDir  \
	--genomeFastaFiles /home/ye/Data/Zoc/Cell/reference/dna/mouse/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa  \
	--sjdbGTFfile /home/ye/Data/Zoc/Cell/reference/dna/mouse/gtf/Mus_musculus.GRCm38.95.gtf \
	--readFilesType Fastx \
	--outSAMtype SAM \
	--runThreadN 2 \
	--outFileNamePrefix /home/ye/Data/Zoc/Cell/reference/data/qc/star/ryanyip_ \
	--readFilesIn /home/ye/Data/Zoc/Cell/reference/data/qc/fastp/B2_L4_A003.R1.clean.fastq.gz   /home/ye/Data/Zoc/Cell/reference/data/qc/fastp/B2_L4_A003.R2.clean.fastq.gz
