import os,re,argparse,datetime
# import ujson as json
import numpy as np
import json

import qtdb_filelist
import extract_func


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extract data from qtdb', allow_abbrev=True)
    parser.add_argument('--path','-p', type=str, help='base path for file read/write (qtdb will be appended)')
    args = parser.parse_args()
    if args.path==None:
        args.path = os.getcwd()


    #get the files in the database
    dbnames=set()
    dbpath = args.path
    targpath = os.path.join(args.path, 'qtdbout')
    if os.path.exists(targpath)==False:
        os.makedirs(targpath)
    db = 'qtdb'
    for f in os.listdir(os.path.join(dbpath,db)):
        fnmatch = re.match('(.*)\.atr$',f)
        if fnmatch:
            fn = fnmatch.group(1)
            dbnames.add( fn )

    all_segments = []
    all_beats = []
    for dbn in dbnames:
        if (dbn in qtdb_filelist.source_details) == False:
            print(dbn + " not found in source detail list")
            continue
        print(dbn)
        all_segments += extract_func.split_qtdb_file(dbpath,dbn,30,qtdb_filelist.source_details[dbn],highpass_freq_cutoff=0.1)
        all_beats += extract_func.pull_qtdb_manual_beats(dbpath,dbn,qtdb_filelist.source_details[dbn],highpass_freq_cutoff=0.1)
        break

    #add our specific classifications to these
    # TODO: check if more than MITDB is required
    # TODO: link the rhythm annotations to the beats
    all_segments = extract_func.classify_mit_segments(all_segments)

    #create the manifest dictionary
    manifest = extract_func.generate_qtdb_manifest(all_segments)
    beat_manifest = extract_func.generate_qtdb_segment_manifest(all_beats)
    with open(os.path.join(targpath,'MANIFEST.json'),'w') as f:
        json.dump(manifest,f,indent=4)
    with open(os.path.join(targpath,'MANIFEST.beats.json'),'w') as f:
        json.dump(beat_manifest,f,indent=4)

    shalist = {}
    for seg in all_segments:
        seg_sha = extract_func.write_atc_from_segment(seg, targpath)
        for k in seg_sha.keys():
            shalist[k]=seg_sha[k]
    for beat in all_beats:
        seg_sha = extract_func.write_beat_atc_from_segment(beat, targpath)
        for k in seg_sha.keys():
            shalist[k]=seg_sha[k]

    with open(os.path.join(targpath,'SHA256SUMS'),'w') as f:
        for s in shalist.keys():
            f.writelines('{}\t{}\n'.format(shalist[s],s))





