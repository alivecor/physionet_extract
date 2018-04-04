import numpy as np
import os,re
import wfdb



def split_file(record_path, record_name, targ_sec):

    #read in the signal and the annotations
    record = wfdb.rdrecord(os.path.join(record_path,record_name))
    recatr = wfdb.rdann(os.path.join(record_path,record_name),'atr')

    if np.any(np.array(record.units)!='mV'):
        raise Exception('alternate units not implemented yet')

    beat_sample = np.array(recatr.sample)
    beat_symbol = np.array(recatr.symbol)

    boundary_ind = np.where(np.array(recatr.symbol)=='+')[0]
    aux_note = [an[1:] for an in np.array(recatr.aux_note)[boundary_ind]]

    #get the boundaries
    boundary_ind = [(boundary_ind[k]+1,boundary_ind[k+1]) for k in range(len(boundary_ind)-1)]
    boundary = [(recatr.sample[b[0]],recatr.sample[b[1]]) for b in boundary_ind]

    atc_span = []
    span_rhythm = []
    for k in range(len(boundary)-1,-1,-1):
        # do we have enough items for a sample?
        if (boundary[k][1]-boundary[k][0])<targ_sec*record.fs:
            boundary.pop(k)
            continue

        num_seg = (boundary[k][1]-boundary[k][0])//(targ_sec*record.fs)
        try_buf = ((boundary[k][1]-boundary[k][0])%(targ_sec*record.fs))//(num_seg+1)

        for br in range(num_seg):
            st = boundary[k][0] + br*(targ_sec*record.fs+try_buf) +try_buf//2
            en = st + targ_sec*record.fs

            #for the current segment, get the midpoint between the first beat within the range and the first beat outside the range
            ix = np.where(beat_sample<st)
            if len(ix[0])!=0 and ix[0][-1]<(len(beat_sample)-1):
                ix = ix[0][-1]
                pre_midpoint = (beat_sample[ix]+beat_sample[ix+1])//2
            else:
                pre_midpoint = st
            ix = np.where(beat_sample>en)
            if len(ix[0])!=0 and ix[0][0]>0:
                ix = ix[0][0]
                post_midpoint = (beat_sample[ix-1]+beat_sample[ix])//2
            else:
                post_midpoint=en

            #average the offsets to try to better center the beats, and if it's less than the try_buf window add it to the ranges
            offset = ((pre_midpoint-st)+(post_midpoint-en))//2
            offset = np.sign(offset)*min(abs(offset),try_buf)
            atc_span.append( (st+offset,en+offset) )
            span_rhythm.append(aux_note[k])


    #now grab all of the beats and sample positions in the sample
    segments = []
    for ind in range(len(atc_span)):
        span=atc_span[ind]
        segments.append({
            'source_ind': span,
            'signal': record.p_signal[span[0]:span[1],:],
            'beat_label': beat_symbol[np.where(np.logical_and(beat_sample>=span[0], beat_sample<span[1]))],
            'beat_index': beat_sample[np.where(np.logical_and(beat_sample>=span[0], beat_sample<span[1]))]-span[0],
            'rhythm_label': span_rhythm[ind],
        })



if __name__ == "__main__":

    #get the files in the database
    dbname=set()
    for f in os.listdir('/Users/schram/projects/physionet_extract/mitdb'):
        fnmatch = re.match('(.*)\.atr$',f)
        if fnmatch:
            fn = fnmatch.group(1)
            dbname.add( 'mitdb/'+fn )


dbname=sorted(list(dbname))
for dbn in dbname:
    recatr = wfdb.rdann('mitdb/'+dbn,'atr')
    print(set(recatr.symbol))


