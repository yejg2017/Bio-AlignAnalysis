import json
import os

def update_json(new_dict,jsonfile=None):
    if jsonfile is None:
        key_dict=new_dict
    else:
        with open(jsonfile,"r") as file:
            key_dict=json.load(file)
            file.close()
        key_dict.update(new_dict)

    jsonfile=jsonfile if jsonfile is not None  else "result.json"

    with open(jsonfile,"w") as file:
        json.dump(key_dict,file)
        file.close()
    print("Update done!")

