import os,re
import ujson as json
import numpy as np

import extract_func

np.seterr(all='raise')


nst_codepath = {
    'db': 'nstdb',
    'split': extract_func.split_physionet_highpass_file,
    'classify': extract_func.classify_mit_segments,
    'manifest': extract_func.generate_nst_manifest,
    'write': extract_func.write_atc_from_segment
}


if __name__ == "__main__":

    code_path = nst_codepath

    #get the files in the database
    dbnames=set()
    dbpath = '/Users/schram/projects/physionet_extract/'
    targpath = '/Users/schram/projects/physionet_extract/nstout'
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
        # all_segments += split_physionet_file(dbpath,db,dbn,30)
        print(dbn)
        all_segments += code_path['split'](dbpath,db,dbn,30)


    #add our specific classifications to these
    all_segments = code_path['classify'](all_segments)

    #create the manifest dictionary
    manifest = code_path['manifest']( all_segments )
    with open(os.path.join(targpath,'MANIFEST.json'),'w') as f:
        json.dump(manifest,f,indent=4)


    shalist = {}
    for seg in all_segments:
        seg_sha = code_path['write'](seg, targpath)
        for k in seg_sha.keys():
            shalist[k]=seg_sha[k]

    with open(os.path.join(targpath,'SHA256SUMS'),'w') as f:
        for s in shalist.keys():
            f.writelines('{}\t{}\n'.format(shalist[s],s))



