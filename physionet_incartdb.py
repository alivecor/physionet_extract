import os,re,argparse
import ujson as json
import numpy as np

import extract_func

np.seterr(all='raise')

import functools
print = functools.partial(print, flush=True)


incartdb_codepath = {
    'db': 'incartdb',
    'split': extract_func.split_physionet_file,
    'classify': extract_func.classify_mit_segments,
    'manifest': extract_func.generate_incartdb_manifest,
    'write': extract_func.write_atc_from_segment
}


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extract data from incartdb', allow_abbrev=True)
    parser.add_argument('--path','-p', type=str, help='base path for file read/write (incartdb will be appended)')
    args = parser.parse_args()
    if args.path==None:
        args.path = os.getcwd()

    code_path = incartdb_codepath

    #get the files in the database
    dbnames=set()
    dbpath = args.path
    targpath = os.path.join(args.path, 'incartdb')
    if os.path.exists(targpath)==False:
        os.makedirs(targpath)
    if os.path.exists(os.path.join(targpath,'data'))==False:
        os.makedirs(os.path.join(targpath,'data'))
    db = code_path['db']
    for f in os.listdir(os.path.join(dbpath,db)):
        fnmatch = re.match('(.*)\.atr$',f)
        if fnmatch:
            fn = fnmatch.group(1)
            dbnames.add( fn )


    all_segments = []
    for dbn in dbnames:
        print(dbn)
        all_segments += code_path['split'](dbpath,db,dbn,30,default_rhythm='N', highpass_freq_cutoff=0.1)


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





