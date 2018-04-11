import os,re
import ujson as json
import numpy as np

import extract_func

np.seterr(all='raise')


cudb_codepath = {
    'db': 'cudb',
    'split': extract_func.split_cudb_file,
    'classify': extract_func.classify_cudb_segments,
    'manifest': extract_func.generate_cudb_manifest,
    'write': extract_func.write_atc_from_segment
}


if __name__ == "__main__":

    code_path = cudb_codepath

    #get the files in the database
    dbnames=set()
    dbpath = '/Users/schram/projects/physionet_extract/'
    targpath = '/Users/schram/projects/physionet_extract/cudbout'
    if os.path.exists(targpath)==False:
        os.makedirs(targpath)
    db = code_path['db']
    for f in os.listdir(os.path.join(dbpath,db)):
        fnmatch = re.match('(.*)\.atr$',f)
        if fnmatch:
            fn = fnmatch.group(1)
            dbnames.add( fn )


    all_segments = []
    for dbn in dbnames:
        if dbn in ('cu12', 'cu15', 'cu24', 'cu25', 'cu32'):
            print(dbn+': PACED; SKIPPING')
            continue
        print(dbn)
        all_segments += code_path['split'](dbpath,db,dbn,30)


    #add our specific classifications to these
    all_segments = code_path['classify'](all_segments)

    #create the manifest dictionary
    manifest = code_path['manifest']( all_segments )
    with open(os.path.join(targpath,'MANIFEST.json'),'w') as f:
        json.dump(manifest,f,indent=4)


    for seg in all_segments:
        code_path['write'](seg, targpath)

