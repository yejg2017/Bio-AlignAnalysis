import os
import sys
import json
import argparse
import numpy as np



def _fastqc(args):
    
    files=[os.path.join(args.fastq_root,file)
            for file in  os.listdir(args.fastq_root) 
            if file.endswith("fastq.gz") or file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fq.gz")]


    for i,file in enumerate(files):
        sh=args.fastqc+" "+"--outdir"+" "+args.outdir
        sh=sh+ " " + "--format"+ " " +args.format
        sh=sh+ " " + "--dir" + " " + args.dir
        sh=sh+ " " + "--threads" + " " + str(args.threads)
        sh=sh+" "+file
        command="nohup "+ sh +"&"
        print("The [%d] commmand is : \n%s"%(i+1,command))
        os.system(command)





if __name__=="__main__":

    parser=argparse.ArgumentParser()
   
    parser.add_argument("--fastq_root",type=str,default=None,help="The root path of fastq files")
    parser.add_argument("--outdir",type=str,default="outdir",help="The output path of fastqc files")
    parser.add_argument("--format",type=str,default="fastq",help="The format of QC file")
    parser.add_argument("--dir",type=str,default="./tmp",help="The temp directory ")
    parser.add_argument("--threads",type=int,default=2,help="The number of threads")
    parser.add_argument("--fastqc",type=str,default="/home/ye/Software/Python/location/envs/Align/bin/fastqc",help="The location path of fastqc") 

    args=parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    
     
    _fastqc(args)
    

