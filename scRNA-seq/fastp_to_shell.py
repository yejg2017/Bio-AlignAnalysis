import os
import sys
import json
import argparse
import numpy as np

#parser=argparse.ArgumentParser()
#parser.add_argument("--in1",type=str,default=None,help="The first fastq input file")
#parser.add_argument("--out1",type=str,default=None,help="The first fastq output file")
#parser.add_argument("--in2",type=str,default=None,help="The second fastq input file")
#parser.add_argument("--out2",type=str,default=None,help="The second fastq output file")
#parser.add_argument("--fastq_root",type=str,default=None,help="The root path of fastq files")
#parser.add_argument("--output",type=str,default="fastp_output",help="The output path of fastp files")
#parser.add_argument("--shell_file",type=str,default=None,help="The shell file to be excuate")

#args=parser.parse_args()


def _read(args):
    
    files=[file for file in  os.listdir(args.fastq_root) 
            if file.endswith("fastq.gz") or \
                    file.endswith(".fastq") or \
                    file.endswith(".fq.gz") or \
                    file.endswith(".fq")]

    unique_groups=[file.split("_")[0]+"_"+file.split("_")[1]  for file in files]
    print("The are :  %d pairefd groups fastq file"%len(np.unique(unique_groups)))

    for group in np.unique(unique_groups):
       file=[file for file in files if group in file]
       assert len(file)==2
       in1,in2=file
       out1=".".join(in1.split(".")[:-2])+".clean.fastq.gz"
       out2=".".join(in2.split(".")[:-2])+".clean.fastq.gz"

       in1=os.path.join(args.fastq_root,in1)
       in2=os.path.join(args.fastq_root,in2)

       out1=os.path.join(args.output,out1)
       out2=os.path.join(args.output,out2)


       shell=args.fastp+" " + "-i" + in1 + \
               " " + "-o" + " " + out1 + \
               " " + "-I" + " " + in2 + \
               " " + "-O" + " " + out2 + \
               " " + "--overrepresentation_analysis"

       print("The  script is : %s \n"%shell)
       os.system("nohup"+" "+shell+"&")



if __name__=="__main__":


    parser=argparse.ArgumentParser()
    #parser.add_argument("--in1",type=str,default=None,help="The first fastq input file")
    #parser.add_argument("--out1",type=str,default=None,help="The first fastq output file")
    #parser.add_argument("--in2",type=str,default=None,help="The second fastq input file")
    #parser.add_argument("--out2",type=str,default=None,help="The second fastq output file")
    parser.add_argument("--fastq_root",type=str,default=None,help="The root path of fastq files")
    parser.add_argument("--output",type=str,default="fastp_output",help="The output path of fastp files")
    parser.add_argument("--fastp",type=str,default="/home/ye/Software/Python/location/envs/Align/bin/fastp",help="The location of fastp")

    args=parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    _read(args)





