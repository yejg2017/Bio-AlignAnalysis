import os
import sys
import json

from .utils import update_json


def _featureCounts(args):
    if args.json is None:
        raise ValueError("Json file can not be None")

    with open(args.json,"r") as file:
        key_dict=json.load(file)
        file.close()

    gtf=key_dict["gtf"]
    samples=list(key_dict["runMode"]["alignReads"]["output"].keys())
 
    for i,sample in enumerate(samples):
        
        if key_dict["sort"]:
            align=key_dict["runMode"]["alignReads"]["sorted_bam"][sample]
        if not key_dict["sort"]:
            align=key_dict["runMode"]["alignReads"]["bam"][sample]

        output=os.path.dirname(align)
        command=args.featureCounts + " " + "-a" + gtf + \
                " " + "-o" + " " + os.path.join(output,"featureCounts.txt") + \
                " " + "-T" + " " + str(args.threads) + \
                " " + "-d" + " " + str(args.min_length) + \
                " " + "-D" + " " + str(args.max_length) + \
                " " + "-g" + " " + "gene_id" +  \
                " " + "-t" + " " + "exon" + \
                " " + "--primary" + \
                " " + align
        log="featureCounts.sh"
        with open(os.path.join(output,log),"w") as f:
           f.write(command)
           f.write("\n")
           f.write("\n")
           whether="if [ $? -ne 0 ];then \n" + \
               "  " + "echo 1 >" + os.path.join(output,"featureCounts.log.out")  + "\n" + \
               "else" + "\n" + \
               "  " + "echo 0 >" + os.path.join(output,"featureCounts.log.out")  + "\n" + \
               "fi"

           cut= "cut -f1,7,8,9,10,11,12"+ " " + os.path.join(output,"featureCounts.txt") + \
                   " " + ">" + " " + os.path.join(output,"featureMatrix.txt")
           f.write(whether)
           f.write("\n")
           f.write("\n")
           f.write(cut)
           f.close()
        os.system("nohup" + " " + "bash" + " " + os.path.join(output,log) + " " + "&")


