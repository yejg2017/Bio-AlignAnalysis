from module import (_STAR,
        genomeGenerate,
        alignReads,
        _featureCounts,
        sam2bam)

import os
import sys
import json
import argparse
import numpy as np


if __name__=="__main__":

    parser=argparse.ArgumentParser("STAR ON EACH GROUP OR ALL,AND THEN FEATURECOUNTS")
    parser.add_argument("--fastq_root",type=str,default="",help="The root path of fastq files")
    
    parser.add_argument("--star",type=str,default="/home/ye/Software/Python/location/envs/Align/bin/STAR",help="The location of STAR")
    parser.add_argument("--runMode",type=str,default="genomeGenerate",help="RunMode : alignReads,genomeGenerate,featureCounts,...")
    parser.add_argument("--genomeDir",type=str,default="/home/ye/Data/Zoc/Cell/reference/data/qc/star/genomeDir",help="The output of genomeDir")
    parser.add_argument("--genomeFastaFiles",type=str,
            default="/home/ye/Data/Zoc/Cell/reference/dna/mouse/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa",
            help="The Fastafiles' path")
    
    parser.add_argument("--sjdbGTFfile",type=str,
            default="/home/ye/Data/Zoc/Cell/reference/dna/mouse/gtf/Mus_musculus.GRCm38.95.gtf",
            help="The path of gtf file")

    parser.add_argument("--readFilesType",type=str,default="Fastx",help="The input file format")
    parser.add_argument("--threads",type=int,default=2,help="The number of threads")
    parser.add_argument("--outTmpKeep",type=str,default="None",help="Whether to keep temp file")
    parser.add_argument("--outSAMtype",type=str,default="SAM",help="The format of output file")
    parser.add_argument("--outFileNamePrefix",type=str,default="/home/ye/Data/Zoc/Cell/reference/data/qc/star")
    

    parser.add_argument("--outFilterMismatchNmax",type=int,default=10,help="The number of outFilterMismatchNma")
    parser.add_argument("--seedSearchStartLmax",type=int,default=50)
    parser.add_argument("--seedPerReadNmax",type=int,default=1000)
    parser.add_argument("--seedPerWindowNmax",type=int,default=50)
    parser.add_argument("--alignTranscriptsPerReadNmax",type=int,default=10000)
    parser.add_argument("--alignTranscriptsPerWindowNmax",type=int,default=100)
    parser.add_argument("--quantMode",type=str,default="GeneCounts")
    
    parser.add_argument("--json",type=str,default=None,help="The result json file")
    
    parser.add_argument("--featureCounts",type=str,default="/home/ye/Software/Python/location/envs/Align/bin/featureCounts",help="where the featureCounts")
    parser.add_argument("--min_length",type=int,default=50)
    parser.add_argument("--max_length",type=int,default=600)

    parser.add_argument("--samtools",type=str,
            default="/home/ye/Software/Python/location/envs/Align/bin/samtools",
            help="where the samtools")
    parser.add_argument("--sort", dest='sort', action='store_true',help="whether sort the bam")
    args=parser.parse_args()

    print("-------------------")
    print("-----STAR Configure-----")
    print("-------------------")
    print("STAR PATH : {}".format(args.star))
    print("runMode : {}".format(args.runMode))
    print("genomeDir : {}".format(args.genomeDir))
    print("genomeFastaFiles : {}".format(args.genomeFastaFiles))
    print("sjdbGTFfile : {}".format(args.sjdbGTFfile))
    print("threads : {}".format(args.threads))
    print("readFilesType : {}".format(args.readFilesType))
    print("outTmpKeep : {}".format(args.outTmpKeep))
    print("outSAMtype : {}".format(args.outSAMtype))
    print("outFileNamePrefix : {}".format(args.outFileNamePrefix))

    print("outFilterMismatchNmax : {}".format(args.outFilterMismatchNmax))
    print("seedSearchStartLmax : {}".format(args.seedSearchStartLmax))
    print("seedPerReadNmax : {}".format(args.seedPerReadNmax))
    print("seedPerWindowNmax : {}".format(args.seedPerWindowNmax))
    print("alignTranscriptsPerReadNmax : {}".format(args.alignTranscriptsPerReadNmax))
    print("alignTranscriptsPerWindowNmax : {}".format(args.alignTranscriptsPerWindowNmax))
    print("quantMode : {}".format(args.quantMode))
    

    print("--------------------")
    print("----featureCounts Configure----")
    print("--------------------")
    print("featureCounts PATH : {}".format(args.featureCounts))
    print("Minimum fragment/template length : {}".format(args.min_length))
    print("Maximum fragment/template length : {}".format(args.max_length))


    print("--------------------")
    print("----samtools Configure----")
    print("--------------------")
    print("samtools PATH : {}".format(args.samtools))
    print("sort : {}".format(args.sort))

    if args.runMode=="alignReads":
    
        alignReads(args)
    
    elif args.runMode=="genomeGenerate":
        genomeGenerate(args)

    elif args.runMode=="samtools":
        print("sam to bam")
        sam2bam(args)
    elif args.runMode=="featureCounts":
        print("featureCounts")
        _featureCounts(args)
    else:
        raise ValueError("This runMode is invalid")

