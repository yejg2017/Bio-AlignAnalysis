import os
import numpy as np
import sys
import json

from .utils  import update_json

def _STAR(args):
    
    if not args.runMode in ["genomeGenerate","alignReads"]:
        raise ValueError("runMode is invalid or try the secondary analysis")

    
    files=[file for file in  os.listdir(args.fastq_root)
            if file.endswith("fastq.gz") or \
                    file.endswith(".fastq") or \
                    file.endswith(".fq.gz") or \
                    file.endswith(".fq")]

    unique_groups=[file.split("_")[0]+"_"+file.split("_")[1]  for file in files]
    print("The are :  %d pairefd groups fastq file"%len(np.unique(unique_groups)))
    
    if args.json is None:
        key_dict={"gtf": args.sjdbGTFfile,
              "fasta": args.genomeFastaFiles,
              "star": args.star,
              "readFilesType": args.readFilesType,
              "runMode" : {
                  "genomeGenerate":{"genomeDir": {},
                                    "log" :      {}},
                  "alignReads"    :{"output": {},
                                    "log" :   {},
                                    "bam" :   {},
                                    "sorted_bam" : {}}
                  },
              "data":{},
              "command" : {"genomeGenerate":{},
                           "alignReads":{}},
              "sort" :args.sort
              }
    else:
        with open(args.json,"r") as file:
            key_dict=json.load(file)
            file.close()
    
    for i,group in enumerate(np.unique(unique_groups)):
       
       file=[file for file in files if group in file]
       assert len(file)==2
       read1=os.path.join(args.fastq_root,file[0])
       read2=os.path.join(args.fastq_root,file[1])
       
       if not os.path.exists(os.path.join(args.outFileNamePrefix,group)):
          os.makedirs(os.path.join(args.outFileNamePrefix,group))

       if not os.path.exists(os.path.join(args.genomeDir,group)):
          os.makedirs(os.path.join(args.genomeDir,group))

       command=args.star + " " + "--runMode" + " " + args.runMode + \
               " " + "--genomeDir" + " " + os.path.join(args.genomeDir,group) + \
               " " + "--genomeFastaFiles" + " " + args.genomeFastaFiles + \
               " " + "--sjdbGTFfile" + " " + args.sjdbGTFfile + \
               " " + "--readFilesType" + " " + args.readFilesType + \
               " " + "--runThreadN" + " " + str(args.threads) + \
               " " + "--outTmpKeep" + " " + args.outTmpKeep + \
               " " + "--outSAMtype" + " " + args.outSAMtype + \
               " " + "--outFileNamePrefix" + " " + os.path.join(args.outFileNamePrefix,group,"raynyip_") + \
               " " + "--outFilterMismatchNmax" + " " + str(args.outFilterMismatchNmax) + \
               " " + "--seedSearchStartLmax" + " " + str(args.seedSearchStartLmax) + \
               " " + "--seedPerReadNmax" + " " + str(args.seedPerReadNmax) + \
               " " + "--seedPerWindowNmax" + " " + str(args.seedPerWindowNmax)+ \
               " " + "--alignTranscriptsPerReadNmax" + " " + str(args.alignTranscriptsPerReadNmax) + \
               " " + "--alignTranscriptsPerWindowNmax" + " " + str(args.alignTranscriptsPerWindowNmax) + \
               " " + "--readFilesCommand" + " " + "zcat" + \
               " " + "--readFilesIn" + " " + read1 + " " + read2
       runMode=args.runMode
       if runMode=="alignReads":
           if args.quantMode is None:
               command=args.star + " " + "--runMode" + " " + args.runMode + \
                   " " + "--genomeDir" + " " + os.path.join(args.genomeDir,group) + \
                   " " + "--sjdbGTFfile" + " " + args.sjdbGTFfile + \
                   " " + "--readFilesType" + " " + args.readFilesType + \
                   " " + "--runThreadN" + " " + str(args.threads) + \
                   " " + "--outTmpKeep" + " " + args.outTmpKeep + \
                   " " + "--outSAMtype" + " " + args.outSAMtype + \
                   " " + "--outFileNamePrefix" + " " + os.path.join(args.outFileNamePrefix,group,"ryanyip_") + \
                   " " + "--outFilterMismatchNmax" + " " + str(args.outFilterMismatchNmax) + \
                   " " + "--seedSearchStartLmax" + " " + str(args.seedSearchStartLmax) + \
                   " " + "--seedPerReadNmax" + " " + str(args.seedPerReadNmax) + \
                   " " + "--seedPerWindowNmax" + " " + str(args.seedPerWindowNmax)+ \
                   " " + "--alignTranscriptsPerReadNmax" + " " + str(args.alignTranscriptsPerReadNmax) + \
                   " " + "--alignTranscriptsPerWindowNmax" + " " + str(args.alignTranscriptsPerWindowNmax) + \
                   " " + "--readFilesCommand" + " " + "zcat" + \
                   " " + "--readFilesIn" + " " + read1 + " " + read2
           else:
                command=args.star + " " + "--runMode" + " " + args.runMode + \
                   " " + "--genomeDir" + " " + os.path.join(args.genomeDir,group) + \
                   " " + "--sjdbGTFfile" + " " + args.sjdbGTFfile + \
                   " " + "--readFilesType" + " " + args.readFilesType + \
                   " " + "--runThreadN" + " " + str(args.threads) + \
                   " " + "--outTmpKeep" + " " + args.outTmpKeep + \
                   " " + "--outSAMtype" + " " + args.outSAMtype + \
                   " " + "--outFileNamePrefix" + " " + os.path.join(args.outFileNamePrefix,group,"ryanyip_") + \
                   " " + "--outFilterMismatchNmax" + " " + str(args.outFilterMismatchNmax) + \
                   " " + "--seedSearchStartLmax" + " " + str(args.seedSearchStartLmax) + \
                   " " + "--seedPerReadNmax" + " " + str(args.seedPerReadNmax) + \
                   " " + "--seedPerWindowNmax" + " " + str(args.seedPerWindowNmax)+ \
                   " " + "--alignTranscriptsPerReadNmax" + " " + str(args.alignTranscriptsPerReadNmax) + \
                   " " + "--alignTranscriptsPerWindowNmax" + " " + str(args.alignTranscriptsPerWindowNmax) + \
                   " " + "--readFilesCommand" + " " + "zcat" + \
                   " " + "--quantMode" + " " + args.quantMode + \
                   " " + "--readFilesIn" + " " + read1 + " " + read2

        
       print("Run Command ")
       #os.system("nohup"+" "+command+"&")

       log=args.runMode+".sh"
       with open(os.path.join(args.outFileNamePrefix,group,log),"w") as f:
           f.write(command)
           f.write("\n")
           f.write("\n")
           whether="if [ $? -ne 0 ];then \n" + \
               "  " + "echo 1 >" + os.path.join(args.outFileNamePrefix,group,args.runMode+".log.out")  + "\n" + \
               "else" + "\n" + \
               "  " + "echo 0 >" + os.path.join(args.outFileNamePrefix,group,args.runMode+".log.out")  + "\n" + \
               "fi"
           f.write(whether)
           f.close()

       #os.system("nohup "+"bash"+ " " + os.path.join(args.outFileNamePrefix,group,log) +"&")
       print("The No. %d  command is %s \n"%(i,command))

       key_dict["data"][group]=[read1,read2]

       if args.runMode=="genomeGenerate":
           key_dict["runMode"]["genomeGenerate"]["genomeDir"][group]=os.path.join(args.genomeDir,group)
           key_dict["runMode"]["genomeGenerate"]["log"][group]=os.path.join(args.outFileNamePrefix,group,args.runMode+".log.out")
           key_dict["command"][args.runMode][group]=os.path.join(args.outFileNamePrefix,group,log)

       if args.runMode=="alignReads":
           key_dict["runMode"]["alignReads"]["output"][group]=os.path.join(args.outFileNamePrefix,group)+"Aligned.out"+"."+args.outSAMtype.lower()
           key_dict["runMode"]["alignReads"]["log"][group]=os.path.join(args.outFileNamePrefix,group,args.runMode+".log.out")
           key_dict["command"][args.runMode][group]=os.path.join(args.outFileNamePrefix,group,log)
    
    update_json(key_dict,args.json)
    #if args.json is None:
    #    jsonfile="result.json"
    #else:
    #    jsonfile=args.json
    #with open(jsonfile,"w") as f:
    #          json.dump(key_dict,f)
    #f.close()




def genomeGenerate(args):
    key_dict={"gtf": args.sjdbGTFfile,
              "fasta": args.genomeFastaFiles,
              "star": args.star,
              "readFilesType": args.readFilesType,
              "runMode" : {
                  "genomeGenerate":{"genomeDir": {},
                                    "log" :      {}},
                  "alignReads"    :{"output": {},
                                    "log" :   {},
                                    "bam" :   {},
                                    "sorted_bam" : {}}
                  },
              "data":{},
              "command" : {"genomeGenerate":{},
                           "alignReads":{}},
              "sort" :args.sort
              }
    if not os.path.exists(args.outFileNamePrefix):
          os.makedirs(args.outFileNamePrefix)

    if not os.path.exists(args.genomeDir):
          os.makedirs(args.genomeDir)

    command=args.star + " " + "--runMode" + " " + "genomeGenerate"  + \
               " " + "--genomeDir" + " " + args.genomeDir + \
               " " + "--genomeFastaFiles" + " " + args.genomeFastaFiles + \
               " " + "--sjdbGTFfile" + " " + args.sjdbGTFfile + \
               " " + "--runThreadN" + " " + str(args.threads) + \
               " " + "--outTmpKeep" + " " + args.outTmpKeep + \
               " " + "--outFileNamePrefix" + " " + os.path.join(args.outFileNamePrefix,"raynyip_") + \
               " " + "--readFilesCommand" + " " + "zcat" 

    print("Command : {}".format(command))
    log="genomeGenerate.sh"
    with open(os.path.join(args.outFileNamePrefix,log),"w") as f:
           f.write(command)
           f.write("\n")
           f.write("\n")
           whether="if [ $? -ne 0 ];then \n" + \
               "  " + "echo 1 >" + os.path.join(args.outFileNamePrefix,"genomeGenerate.log.out")  + "\n" + \
               "else" + "\n" + \
               "  " + "echo 0 >" + os.path.join(args.outFileNamePrefix,"genomeGenerate.log.out")  + "\n" + \
               "fi"
           f.write(whether)
           f.close()
    key_dict["runMode"]["genomeGenerate"]["genomeDir"]=args.genomeDir
    key_dict["runMode"]["genomeGenerate"]["log"]=os.path.join(args.outFileNamePrefix,"genomeGenerate.log.out")
    key_dict["command"]["genomeGenerate"]=os.path.join(args.outFileNamePrefix,log)

    os.system("nohup "+"bash"+ " " + os.path.join(args.outFileNamePrefix,log) +"&")
    update_json(key_dict,args.json)

def alignReads(args):
    
    if args.json is None:
        key_dict={"gtf": args.sjdbGTFfile,
              "fasta": args.genomeFastaFiles,
              "star": args.star,
              "readFilesType": args.readFilesType,
              "runMode" : {
                  "genomeGenerate":{"genomeDir": {},
                                    "log" :      {}},
                  "alignReads"    :{"output": {},
                                    "log" :   {},
                                    "bam" :   {},
                                    "sorted_bam" : {}}
                  },
              "data":{},
              "command" : {"genomeGenerate":{},
                           "alignReads":{}},
              "sort" :args.sort
              }
    else:
        with open(args.json,"r") as file:
            key_dict=json.load(file)
            file.close()

    
    files=[file for file in  os.listdir(args.fastq_root)
            if file.endswith("fastq.gz") or \
                    file.endswith(".fastq") or \
                    file.endswith(".fq.gz") or \
                    file.endswith(".fq")]
    unique_groups=[file.split("_")[0]+"_"+file.split("_")[1]  for file in files]
    print("The are :  %d pairefd groups fastq file"%len(np.unique(unique_groups)))


    for i,group in enumerate(np.unique(unique_groups)):

       file=[file for file in files if group in file]
       assert len(file)==2
       read1=os.path.join(args.fastq_root,file[0])
       read2=os.path.join(args.fastq_root,file[1])

       if not os.path.exists(os.path.join(args.outFileNamePrefix,group)):
          os.makedirs(os.path.join(args.outFileNamePrefix,group))

       if not os.path.exists(os.path.join(args.genomeDir,group)):
          os.makedirs(os.path.join(args.genomeDir,group))

       
       if args.quantMode is None:
               command=args.star + " " + "--runMode" + " " + "alignReads" + \
                   " " + "--genomeDir" + " " + key_dict["runMode"]["genomeGenerate"]["genomeDir"] + \
                   " " + "--sjdbGTFfile" + " " + args.sjdbGTFfile + \
                   " " + "--readFilesType" + " " + args.readFilesType + \
                   " " + "--runThreadN" + " " + str(args.threads) + \
                   " " + "--outTmpKeep" + " " + args.outTmpKeep + \
                   " " + "--outSAMtype" + " " + args.outSAMtype + \
                   " " + "--outFileNamePrefix" + " " + os.path.join(args.outFileNamePrefix,group,"ryanyip_") + \
                   " " + "--outFilterMismatchNmax" + " " + str(args.outFilterMismatchNmax) + \
                   " " + "--seedSearchStartLmax" + " " + str(args.seedSearchStartLmax) + \
                   " " + "--seedPerReadNmax" + " " + str(args.seedPerReadNmax) + \
                   " " + "--seedPerWindowNmax" + " " + str(args.seedPerWindowNmax)+ \
                   " " + "--alignTranscriptsPerReadNmax" + " " + str(args.alignTranscriptsPerReadNmax) + \
                   " " + "--alignTranscriptsPerWindowNmax" + " " + str(args.alignTranscriptsPerWindowNmax) + \
                   " " + "--genomeLoad" + " " + "LoadAndKeep" + \
                   " " + "--readFilesCommand" + " " + "zcat" + \
                   " " + "--readFilesIn" + " " + read1 + " " + read2
       else:
                command=args.star + " " + "--runMode" + " " + args.runMode + \
                   " " + "--genomeDir" + " " + key_dict["runMode"]["genomeGenerate"]["genomeDir"] + \
                   " " + "--sjdbGTFfile" + " " + args.sjdbGTFfile + \
                   " " + "--readFilesType" + " " + args.readFilesType + \
                   " " + "--runThreadN" + " " + str(args.threads) + \
                   " " + "--outTmpKeep" + " " + args.outTmpKeep + \
                   " " + "--outSAMtype" + " " + args.outSAMtype + \
                   " " + "--outFileNamePrefix" + " " + os.path.join(args.outFileNamePrefix,group,"ryanyip_") + \
                   " " + "--outFilterMismatchNmax" + " " + str(args.outFilterMismatchNmax) + \
                   " " + "--seedSearchStartLmax" + " " + str(args.seedSearchStartLmax) + \
                   " " + "--seedPerReadNmax" + " " + str(args.seedPerReadNmax) + \
                   " " + "--seedPerWindowNmax" + " " + str(args.seedPerWindowNmax)+ \
                   " " + "--alignTranscriptsPerReadNmax" + " " + str(args.alignTranscriptsPerReadNmax) + \
                   " " + "--alignTranscriptsPerWindowNmax" + " " + str(args.alignTranscriptsPerWindowNmax) + \
                   " " + "--genomeLoad" + " " + "LoadAndKeep" + \
                   " " + "--readFilesCommand" + " " + "zcat" + \
                   " " + "--quantMode" + " " + args.quantMode + \
                   " " + "--readFilesIn" + " " + read1 + " " + read2
       log="alignReads.sh"
       with open(os.path.join(args.outFileNamePrefix,group,log),"w") as f:
            f.write(command)
            f.write("\n")
            f.write("\n")
            whether="if [ $? -ne 0 ];then \n" + \
               "  " + "echo 1 >" + os.path.join(args.outFileNamePrefix,group,"alignReads.log.out")  + "\n" + \
               "else" + "\n" + \
               "  " + "echo 0 >" + os.path.join(args.outFileNamePrefix,group,"alignReads.log.out")  + "\n" + \
               "fi"
            f.write(whether)
            f.close()

       os.system("nohup "+"bash"+ " " + os.path.join(args.outFileNamePrefix,group,log) +"&")
       print("The No. %d  command is %s \n"%(i,command))
       print("Update json file...")
       key_dict["data"][group]=[read1,read2]

       key_dict["runMode"]["alignReads"]["output"][group]=os.path.join(args.outFileNamePrefix,group,"Aligned.out"+"."+args.outSAMtype.lower())
       key_dict["runMode"]["alignReads"]["log"][group]=os.path.join(args.outFileNamePrefix,group,"alignReads.log.out")
       key_dict["command"][args.runMode][group]=os.path.join(args.outFileNamePrefix,group,log)

    update_json(key_dict,args.json)



