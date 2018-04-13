import numpy as np
import os,re,hashlib
import ujson as json
from scipy import signal
import scipy.signal
import wfdb
from acpy.atc import write_atc_file


atcfs=300

def gen_atc_filebase(seg,ind_lead):
    return '{}_{}_{}'.format(seg['source_file'],seg['source_ind'][0],seg['header_leads'][ind_lead])

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
    boundary_ind = [(raw_bound[k],raw_bound[k+1]) for k in range(len(raw_bound)-1)]
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

            atc_span.append( (st,en) )
            span_rhythm.append(aux_note[k])

    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]

        #rescale and resample the actual signal
        numerical_error = False
        rescaled_signal = np.zeros(shape=[record.p_signal.shape[1],atcfs*split_length], dtype=np.int16)
        for il in range(record.p_signal.shape[1]):
            if np.any(np.isnan(record.p_signal[span[0]:span[1],il])):
                print('nan found in {} at sample index {}; skipping'.format(record_name,span))
                numerical_error = True
                continue

            resamp = 2000*signal.resample_poly(record.p_signal[span[0]:span[1],il],atcfs,record.fs)

            if np.any(np.abs(resamp)>=32768):
                print('int16 overflow in {} in sample range {}; skipping'.format(record_name,span))
                numerical_error = True
                continue
                # raise Exception('int16 overflow')
            rescaled_signal[il,:] = resamp.astype(np.int16)
        if numerical_error:
            continue

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


#split_cudb_file - this function has similar logic to the split_physionet_file code, but has additional logic to assume that the 
# initial rhythm for each file is normal (rhythm annotations are more sparse in these files) 
def split_cudb_file(record_path, dbname, record_name, split_length, verbose=True):
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
    aux_note.insert(0, 'N')         #assume the rhythm is initially normal

    if raw_bound.size > 0:
        #get the boundaries
        boundary_ind = [(raw_bound[k],raw_bound[k+1]) for k in range(len(raw_bound)-1)]
        boundary_ind.insert( 0, (0, raw_bound[0]) )
        boundary_ind.append( (raw_bound[len(raw_bound)-1],-1) )
        boundary = [(recatr.sample[b[0]],recatr.sample[b[1]]) for b in boundary_ind]
    else:
        boundary = [(0,record.p_signal.shape[0])]

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

            #make sure that we have beats and that the gap doesn't exceed 2s (indicative of unlabeled arrhythmia)
            beats = beat_sample[np.where( np.logical_and(beat_sample>=st, beat_sample<en) )]
            if beats.size==0:
                print('no beats in {}: {}'.format(record_name, (st,en)))
                continue
            if beats[0]>(st+2*record.fs) or beats[-1]<(en-2*record.fs) or np.any( (beats[1:]-beats[:-1])>2*record.fs ):
                print('excessively long gap between beats in {}: {}'.format(record_name, (st,en)))
                continue

            atc_span.append( (st,en) )
            span_rhythm.append( aux_note[k] )

    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]

        #rescale and resample the actual signal
        overflow = False
        rescaled_signal = np.zeros(shape=[record.p_signal.shape[1],atcfs*split_length], dtype=np.int16)
        for il in range(record.p_signal.shape[1]):
            if np.any(np.isnan(record.p_signal[span[0]:span[1],il])):
                print('nan found in {} at sample index {}; skipping'.format(record_name,span))
                numerical_error = True
                continue

            resamp = 2000*signal.resample_poly(record.p_signal[span[0]:span[1],il],atcfs,record.fs)

            if np.any(np.abs(resamp)>=32768):
                print('int16 overflow in {} at sample index {}; skipping'.format(record_name,span[0]))
                overflow = True
                continue
                # raise Exception('int16 overflow')
            rescaled_signal[il,:] = resamp.astype(np.int16)
        if overflow:
            continue

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


#split_vfdb_file - this function has similar logic to the split_physionet_file code, but does not include any logic involving heartbeats
# as there are no beat annotations provided in vfdb
def split_vfdb_file(record_path, dbname, record_name, split_length, verbose=True):
    #read in the signal and the annotations
    record = wfdb.rdrecord(os.path.join(record_path,dbname,record_name))
    recatr = wfdb.rdann(os.path.join(record_path,dbname,record_name),'atr')

    #read in the header and get the leads
    with open(os.path.join(record_path,dbname,record_name+'.hea'),'r') as headf:
        headf.readline()
        header_leads = [headf.readline().strip().split(' ')[-1]+str(il) for il in range(record.p_signal.shape[1]) ]


    if np.any(np.array(record.units)!='mV'):
        raise Exception('alternate units not implemented yet')

    beat_symbol = np.array(recatr.symbol)

    raw_bound = np.where(np.array(recatr.symbol)=='+')[0]
    aux_note = [  an[1:].split('\x00')[0] for an in np.array(recatr.aux_note)[raw_bound]]

    #get the boundaries
    boundary_ind = [(raw_bound[k],raw_bound[k+1]) for k in range(len(raw_bound)-1)]
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

            atc_span.append( (st,en) )
            span_rhythm.append( aux_note[k] )

    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]

        #rescale and resample the actual signal
        overflow = False
        rescaled_signal = np.zeros(shape=[record.p_signal.shape[1],atcfs*split_length], dtype=np.int16)
        for il in range(record.p_signal.shape[1]):
            if np.any(np.isnan(record.p_signal[span[0]:span[1],il])):
                print('nan found in {} at sample index {}; skipping'.format(record_name,span))
                numerical_error = True
                continue

            resamp = 2000*signal.resample_poly(record.p_signal[span[0]:span[1],il],atcfs,record.fs)

            if np.any(np.abs(resamp)>=32768):
                print('int16 overflow in {} at sample index {}; skipping'.format(record_name,span[0]))
                overflow = True
                continue
                # raise Exception('int16 overflow')
            rescaled_signal[il,:] = resamp.astype(np.int16)
        if overflow:
            continue

        segments.append({
            'source_file': dbname+record_name,
            'source_ind': span,
            'signal': rescaled_signal,
            'beat_label': [],
            'beat_index': [],
            'rhythm_label': span_rhythm[ind],
            'header_leads': header_leads
        })

    sorted(segments, key=lambda x: x['source_ind'][0])

    return(segments)


#split_afdb_file - this function has similar logic to the split_physionet_file code, but changes how the 
# boundaries are calculated as there are only rhythm annotations in the atr files in afdb, and the beat annotations are in the qrs files
def split_afdb_file(record_path, dbname, record_name, split_length):
    #read in the signal and the annotations
    record = wfdb.rdrecord(os.path.join(record_path,dbname,record_name))
    recatr = wfdb.rdann(os.path.join(record_path,dbname,record_name),'atr')
    recqrs = wfdb.rdann(os.path.join(record_path,dbname,record_name),'qrs')

    #read in the header and get the leads
    with open(os.path.join(record_path,dbname,record_name+'.hea'),'r') as headf:
        headf.readline()
        header_leads = [headf.readline().strip().split(' ')[-1] for il in range(record.p_signal.shape[1]) ]


    if np.any(np.array(record.units)!='mV'):
        raise Exception('alternate units not implemented yet')

    beat_sample = np.array(recqrs.sample)
    beat_symbol = np.array(recqrs.symbol)

    raw_bound = np.where(np.array(recatr.symbol)=='+')[0]
    aux_note = [  an[1:].split('\x00')[0] for an in np.array(recatr.aux_note)[raw_bound]]

    #get the boundaries
    boundary_ind = [(raw_bound[k],raw_bound[k+1]) for k in range(len(raw_bound)-1)]
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

            atc_span.append( (st,en) )
            span_rhythm.append(aux_note[k])

    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]

        #rescale and resample the actual signal
        numerical_error = False
        rescaled_signal = np.zeros(shape=[record.p_signal.shape[1],atcfs*split_length], dtype=np.int16)
        for il in range(record.p_signal.shape[1]):
            if np.any(np.isnan(record.p_signal[span[0]:span[1],il])):
                print('nan found in {} at sample index {}; skipping'.format(record_name,span))
                numerical_error = True
                continue

            resamp = 2000*signal.resample_poly(record.p_signal[span[0]:span[1],il],atcfs,record.fs)

            if np.any(np.abs(resamp)>=32768):
                print('int16 overflow in {} in sample range {}; skipping'.format(record_name,span))
                numerical_error = True
                continue
                # raise Exception('int16 overflow')
            rescaled_signal[il,:] = resamp.astype(np.int16)
        if numerical_error:
            continue

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


#the same logic as to split_physionet_filter, with a 1st order butterworth highpass filter applied to the signal (applied twice, per scipy.signal.filtfilt)
def split_physionet_highpass_file(record_path, dbname, record_name, split_length, freq_cutoff=0.1):
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
    boundary_ind = [(raw_bound[k],raw_bound[k+1]) for k in range(len(raw_bound)-1)]
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

    # apply a 0.1Hz highpass butterworth filter to the signal
    b,a=scipy.signal.butter(N=1,Wn=freq_cutoff*np.pi/record.fs,btype='highpass',analog=False)
    filt_signal = scipy.signal.filtfilt(b,a, record.p_signal, axis=0)

    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]

        #rescale and resample the actual signal
        overflow = False
        rescaled_signal = np.zeros(shape=[filt_signal.shape[1],atcfs*split_length], dtype=np.int16)
        for il in range(filt_signal.shape[1]):
            resamp = 2000*signal.resample_poly(filt_signal[span[0]:span[1],il],atcfs,record.fs)

            if np.any(np.abs(resamp)>=32768):
                print('int16 overflow in {} at sample index {}; skipping'.format(record_name,span[0]))
                overflow = True
                continue
                # raise Exception('int16 overflow')
            rescaled_signal[il,:] = resamp.astype(np.int16)
        if overflow:
            continue

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


def write_atc_from_segment(seg, targpath):

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

    shalist = {}
    for il in range(seg['signal'].shape[0]):
        atc_data_dict['ecg']={}
        atc_data_dict['ecg']['DataLen'] = 2*seg['signal'].shape[1]
        atc_data_dict['ecg']['Data'] = seg['signal'][il,:]

        atcfn = gen_atc_filebase(seg,il)+'.atc'
        write_atc_file.write_atc_file(os.path.join(targpath,atcfn), atc_data_dict)

        #calculate the sha256 of the atc just written
        hf = hashlib.sha256()
        with open(os.path.join(targpath,atcfn),'rb') as f:
            hf.update( f.read() )
        shalist[atcfn] = hf.hexdigest()

        #write an additonal json with the algsuite_target annotation
        algsuite_target = {
            'globalLabels': seg['alg_label']
        }
        algjsonfn = gen_atc_filebase(seg,il)+'.algsuite_target.json'
        with open(os.path.join(targpath,algjsonfn),'w') as f:
            json.dump(algsuite_target,f)

        #calculate the sha256 of the atc just written
        hf = hashlib.sha256()
        with open(os.path.join(targpath,algjsonfn),'rb') as f:
            hf.update( f.read() )
        shalist[algjsonfn] = hf.hexdigest()

    return shalist



def add_segments_to_manifest(manifest, seg_list):
    for seg in seg_list:
        for il in range(len(seg['header_leads'])):
            manifest['recordings'].append({
                'filename': gen_atc_filebase(seg,il)+'.atc',
                'type': 'atc',
                'metadata': {
                    'lead': seg['header_leads'][il],
                    'source_record': seg['source_file'],
                    'source_index': seg['source_ind']
                }
            })

    return manifest


def generate_mit_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'mitdb',
            'description': 'MIT-BIH Arrhythmia Database',
            'source': 'Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001). (PMID: 11446209)'
        },
        'recordings': [],
        'labels': {
            'algsuite_target':{
                'description': 'the result we expect alg-suite to produce',
                'type': 'rhythm'
            }
        }
    }
    add_segments_to_manifest(manifest, seg_list)
    return manifest


def generate_nst_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'nstdb',
            'description': 'MIT-BIH Noise Stress Test Database, sent through 0.1Hz butterworth filter. Excluded samples with signal above atc range (+/-16mV)',
            'source': 'Moody GB, Muldrow WE, Mark RG. A noise stress test for arrhythmia detectors. Computers in Cardiology 1984; 11:381-384'
        },
        'recordings': [],
        'labels': {
            'algsuite_target':{
                'description': 'the result we expect alg-suite to produce',
                'type': 'rhythm'
            }
        }
    }
    add_segments_to_manifest(manifest, seg_list)
    return manifest


def generate_edb_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'edb',
            'description': 'European ST-T Database',
            'source': 'Taddei A, Distante G, Emdin M, Pisani P, Moody GB, Zeelenberg C, Marchesi C. The European ST-T Database: standard for evaluating systems for the analysis of ST-T changes in ambulatory electrocardiography. European Heart Journal 13: 1164-1172 (1992)'
        },
        'recordings': [],
        'labels': {
            'algsuite_target':{
                'description': 'the result we expect alg-suite to produce',
                'type': 'rhythm'
            }
        }
    }
    add_segments_to_manifest(manifest, seg_list)
    return manifest


def generate_cudb_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'cudb',
            'description': 'Creighton University Ventricular Tachyarrhythmia Database',
            'source': 'Nolle FM, Badura FK, Catlett JM, Bowser RW, Sketch MH. CREI-GARD, a new concept in computerized arrhythmia monitoring systems. Computers in Cardiology 13:515-518 (1986)',
            'notes': 'Matt is suspicious of the rhythm labels or lack thereof'
        },
        'recordings': [],
        'labels': {
            'algsuite_target':{
                'description': 'the result we expect alg-suite to produce',
                'type': 'rhythm'
            }
        }
    }
    add_segments_to_manifest(manifest, seg_list)
    return manifest


def generate_afdb_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'afdb',
            'description': 'MIT-BIH Atrial Fibrillation Database',
            'source': 'Moody GB, Mark RG. A new method for detecting atrial fibrillation using R-R intervals. Computers in Cardiology. 10:227-230 (1983)',
        },
        'recordings': [],
        'labels': {
            'algsuite_target':{
                'description': 'the result we expect alg-suite to produce',
                'type': 'rhythm'
            }
        }
    }
    add_segments_to_manifest(manifest, seg_list)
    return manifest


def generate_vfdb_manifest(seg_list):
    manifest = {
        'dataset': {
            'name': 'vfdb',
            'description': 'MIT-BIH Malignant Ventricular Arrhythmia Database',
            'source': 'Greenwald SD. Development and analysis of a ventricular fibrillation detector. M.S. thesis, MIT Dept. of Electrical Engineering and Computer Science, 1986',
            'notes': 'No beat labels in this dataset'
        },
        'recordings': [],
        'labels': {
            'algsuite_target':{
                'description': 'the result we expect alg-suite to produce',
                'type': 'rhythm'
            }
        }
    }
    add_segments_to_manifest(manifest, seg_list)
    return manifest


def classify_mit_segments(all_segments):
    # filter out any segments that are paced
    all_segments = [s for s in all_segments if s['rhythm_label']!='P']

    #apply the desired alg-suite label for each segment
    for seg in all_segments:
        if seg['rhythm_label'] == 'VT' or seg['rhythm_label'] == 'VFL':
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


def classify_cudb_segments(all_segments):
    # filter out any segments that are paced
    all_segments = [s for s in all_segments if s['rhythm_label']!='P']

    #apply the desired alg-suite label for each segment
    for seg in all_segments:
        if seg['rhythm_label'] == 'VT' or seg['rhythm_label'] == 'VF':
            seg['alg_label'] = ['unclassified']
            continue
        if seg['rhythm_label'] == 'AF':
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



def classify_vfdb_segments(all_segments):
    # filter out any segments that are paced
    all_segments = [s for s in all_segments if s['rhythm_label']!='P']

    #apply the desired alg-suite label for each segment
    for seg in all_segments:
        #ventricular fibrillation, tachycardia, or fibrillation
        if seg['rhythm_label'] == 'VT' or seg['rhythm_label'] == 'VF' or seg['rhythm_label'] == 'VFIB' or seg['rhythm_label'] == 'VFL':
            seg['alg_label'] = ['unclassified']
            continue

        #normal or unclassified?
        if seg['rhythm_label'] == 'N' or seg['rhythm_label']=='SBR' or seg['rhythm_label']=='NSR':
            seg['alg_label'] = ['normal']
        if seg['rhythm_label'] == 'NOISE':
            seg['alg_label'] = ['noise']
        else:
            seg['alg_label'] = ['unclassified']


    return(all_segments)

