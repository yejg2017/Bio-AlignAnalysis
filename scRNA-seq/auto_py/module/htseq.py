import os
import sys
import json

from .utils import update_json


def _htseq_count(args):
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
        command=args.htseq+ " " + "--format" + " " + "bam" + \
                " " + "-s" + " " + "no" + \
                " " + "-t" + " " + "exon" + \
                " " + "-i" + " " + "gene_id" + \
                " " + "--additional-attr=gene_name"  + \
                " " + align + " " + args.gff + " " + ">" + " " + os.path.join(output,"htseq-count.txt")

        log="HTSeqCounts.sh"
        with open(os.path.join(output,log),"w") as f:
           f.write(command)
           f.write("\n")
           f.write("\n")
           whether="if [ $? -ne 0 ];then \n" + \
               "  " + "echo 1 >" + os.path.join(output,"htseq-count.log.out")  + "\n" + \
               "else" + "\n" + \
               "  " + "echo 0 >" + os.path.join(output,"htseq-count.log.out")  + "\n" + \
               "fi"
           f.write(whether)
           f.close()

        print("The No. %d  command is %s \n"%(i,command))
        os.system("nohup" + " " + "bash" + " " + os.path.join(output,log) + " " + "&")
        #update_json(args.json,key_dict)

