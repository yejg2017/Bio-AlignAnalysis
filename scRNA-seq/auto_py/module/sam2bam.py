import os
import sys
import json

from .utils import update_json


           
def sam2bam(args):
    if args.json is None:
        raise ValueError("Json file can not be None")

    with open(args.json,"r") as file:
        key_dict=json.load(file)
        file.close()
    

    samples=list(key_dict["runMode"]["alignReads"]["output"].keys())
    for i,sample in enumerate(samples):

        align=key_dict["runMode"]["alignReads"]["output"][sample]
        output=os.path.join(os.path.dirname(align),sample)

        command=args.samtools + " " + "view" + \
                " " + "-S" + " " + \
                " " + "-b" + " " + \
                " " + align + " " +  \
                " " + ">" + " " + \
                " " + os.path.join(output,"Align.out.bam")
        key_dict["runMode"]["alignReads"]["bam"][sample]=os.path.join(output,"Align.out.bam")

        if args.sort:
            bam_sort=args.samtools + " " + "sort" + \
                    " " + os.path.join(output,"Align.out.bam") + \
                    " " + "-o" + \
                    " " + os.path.join(output,"Align.out.sorted.bam")
            key_dict["runMode"]["alignReads"]["sorted_bam"][sample]=os.path.join(output,"Align.out.sorted.bam")
        else:
            bam_sort=""

        log="sam2bam.sh"
        with open(os.path.join(output,log),"w") as f:
           f.write(command)
           f.write("\n")
           f.write(bam_sort)
           f.write("\n")
           f.write("\n")
           whether="if [ $? -ne 0 ];then \n" + \
               "  " + "echo 1 >" + os.path.join(output,"sam2bam.log.out")  + "\n" + \
               "else" + "\n" + \
               "  " + "echo 0 >" + os.path.join(output,"sam2bam.log.out")  + "\n" + \
               "fi"
           f.write(whether)
           f.close()

        print("The No. %d  command is %s \n"%(i,command))
        os.system("nohup" + " " + "bash" + " " + os.path.join(output,log) + " " + "&")
        update_json(key_dict,args.json)
          
