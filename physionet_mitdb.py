import numpy as np
import os,re
import ujson as json
from scipy import signal
import wfdb
from acpy.atc import write_atc_file

np.seterr(all='raise')


atcfs=300

def gen_atc_filename(seg,ind_lead):
    return '{}_{}_{}.atc'.format(seg['source_file'],seg['source_ind'][0],seg['header_leads'][ind_lead])

def split_physionet_file(record_path, dbname, record_name, split_length):
    #read in the signal and the annotations
    record = wfdb.rdrecord(os.path.join(record_path,dbname,record_name))
    recatr = wfdb.rdann(os.path.join(record_path,dbname,record_name),'atr')

    #read in the header and get the leads
    with open(os.path.join(record_path,dbname,record_name+'.hea'),'r') as headf:
        headf.readline()
        header_leads = [headf.readline().strip().split(' ')[-1] for il in range(record.p_signal.shape[1]) ]


    if np.any(np.array(record.units)!='mV'):
        raise Exception('alternate units not implemented yet')

    beat_sample = np.array(recatr.sample)
    beat_symbol = np.array(recatr.symbol)

    raw_bound = np.where(np.array(recatr.symbol)=='+')[0]
    aux_note = [  an[1:].split('\x00')[0] for an in np.array(recatr.aux_note)[raw_bound]]

    #get the boundaries
    boundary_ind = [(raw_bound[k]+1,raw_bound[k+1]) for k in range(len(raw_bound)-1)]
    boundary_ind.append( (raw_bound[len(raw_bound)-1],-1) )
    boundary = [(recatr.sample[b[0]],recatr.sample[b[1]]) for b in boundary_ind]

    atc_span = []
    span_rhythm = []
    for k in range(len(boundary)-1,-1,-1):
        # do we have enough items for a sample?
        if (boundary[k][1]-boundary[k][0])<split_length*record.fs:
            boundary.pop(k)
            continue

        num_seg = (boundary[k][1]-boundary[k][0])//(split_length*record.fs)

        for br in range(num_seg):
            st = int(boundary[k][0] + br*(split_length*record.fs))
            en = int(st + split_length*record.fs)

            #average the offsets to try to better center the beats, and if it's less than the try_buf window add it to the ranges
            atc_span.append( (st,en) )
            span_rhythm.append(aux_note[k])


    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]

        #rescale and resample the actual signal
        rescaled_signal = np.zeros(shape=[record.p_signal.shape[1],atcfs*split_length], dtype=np.int16)
        for il in range(record.p_signal.shape[1]):
            resamp = 2000*signal.resample_poly(record.p_signal[span[0]:span[1],il],atcfs,record.fs)
            if np.any(np.abs(resamp)>=32768):
                raise Exception('int16 overflow')
            rescaled_signal[il,:] = resamp.astype(np.int16)

        #rescale the beat indices
        beat_labels = beat_symbol[np.where(np.logical_and(beat_sample>=span[0], beat_sample<span[1]))]
        beat_samples = ((beat_sample[np.where(np.logical_and(beat_sample>=span[0], beat_sample<span[1]))]-span[0])*atcfs/record.fs).astype(np.uint32)

        segments.append({
            'source_file': dbname+record_name,
            'source_ind': span,
            'signal': rescaled_signal,
            'beat_label': beat_labels,
            'beat_index': beat_samples,
            'rhythm_label': span_rhythm[ind],
            'header_leads': header_leads
        })

    sorted(segments, key=lambda x: x['source_ind'][0])

    return(segments)

def write_atc_from_segment(seg, dbpath):

    atc_data_dict={}

    # File signature and file version
    atc_data_dict['header'] = {}
    atc_data_dict['header']['FileSig'] = (b'ALIVE', 0, 0, 0)
    atc_data_dict['header']['FileVer'] = 2

    # Info Block
    atc_data_dict['info'] = {}
    atc_data_dict['info']['DataLen'] = 264
    atc_data_dict['info']['DateRec'] = '1975-01-01T00:00:00.000+00:00'.encode('utf-8')
    atc_data_dict['info']['RecUUID'] = ''.encode('utf-8')
    atc_data_dict['info']['PhoneUDID'] = ''.encode('utf-8')
    atc_data_dict['info']['PhoneModel'] = ''.encode('utf-8')
    atc_data_dict['info']['RecSW'] = ''.encode('utf-8')
    atc_data_dict['info']['RecHW'] = 'holter'.encode('utf-8')
    atc_data_dict['info']['Loc'] = ''.encode('utf-8')

    # ECG Format Block
    atc_data_dict['fmt'] = {}
    atc_data_dict['fmt']['DataLen'] = 8
    atc_data_dict['fmt']['ECGFormat'] = 1
    atc_data_dict['fmt']['Fs'] = 300
    atc_data_dict['fmt']['AmpRes_nV'] = 500
    atc_data_dict['fmt']['Flags'] = 0
    atc_data_dict['fmt']['Reserved'] = 0

    for il in range(seg['signal'].shape[0]):
        atc_data_dict['ecg']={}
        atc_data_dict['ecg']['DataLen'] = 2*seg['signal'].shape[1]
        atc_data_dict['ecg']['Data'] = seg['signal'][il,:]

        write_atc_file.write_atc_file(os.path.join(dbpath,gen_atc_filename(seg,il)), atc_data_dict)


def generate_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'mitdb',
            'description': 'MIT-BIH Arrhythmia Database',
            'source': 'Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001). (PMID: 11446209)'
        },
        'recordings': []
    }
    for seg in seg_list:
        for il in range(len(seg['header_leads'])):
            manifest['recordings'].append({
                'filename': gen_atc_filename(seg,il),
                'type': 'atc',
                'metadata': {
                    'lead': seg['header_leads'][il],
                    'source_record': seg['source_file'],
                    'source_index': seg['source_ind'],
                    'algsuite_target': seg['alg_label']
                }
            })

    return manifest


def classify_mit_segments(all_segments):
    # filter out any segments that are paced
    all_segments = [s for s in all_segments if s['rhythm_label']!='P']

    #apply the desired alg-suite label for each segment
    for seg in all_segments:
        if seg['rhythm_label'] == 'VT':
            seg['alg_label'] = ['unclassified']
            continue
        if seg['rhythm_label'] == 'AFIB' or seg['rhythm_label'] == 'AFL':
            seg['alg_label'] = ['afib']
            continue

        #compute the heart rate from the median difference between beats
        median_interval = np.median(seg['beat_index'][1:]-seg['beat_index'][:-1])
        heart_rate = 60*atcfs/median_interval

        #are we outside our acceptable heart beat range?
        if heart_rate<40 or heart_rate>140:
            seg['alg_label'] = ['unclassified']
            continue

        #normal or unclassified?
        if seg['rhythm_label'] == 'N' or seg['rhythm_label']=='SBR':
            seg['alg_label'] = ['normal']
        else:
            seg['alg_label'] = ['unclassified']

        #classify brady or tachycardia
        if heart_rate<50:
            seg['alg_label'].append('brady')
            continue
        if heart_rate>=100:
            seg['alg_label'].append('tachy')
            continue

    return(all_segments)


mit_codepath = {
    'db': 'mitdb',
    'split': split_physionet_file,
    'classify': classify_mit_segments,
    'manifest': generate_manifest,
    'write': write_atc_from_segment
}

if __name__ == "__main__":

    code_path = mit_codepath

    #get the files in the database
    dbnames=set()
    dbpath = '/Users/schram/projects/physionet_extract/'
    targpath = '/Users/schram/projects/physionet_extract/mitout'
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
        all_segments += code_path['split'](dbpath,db,dbn,30)


    #add our specific classifications to these
    all_segments = code_path['classify'](all_segments)

    #create the manifest dictionary
    manifest = code_path['manifest']( all_segments )
    with open(os.path.join(targpath,'MANIFEST.json'),'w') as f:
        json.dump(manifest,f,indent=4)


    for seg in all_segments:
        code_path['write'](seg, targpath)





