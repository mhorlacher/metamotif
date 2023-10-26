# %%
import copy

import gin
import numpy as np

from metamotif.utils import running_mean

# %%
def find_significance_threshold(scores, size=2, n=1000, sig_p=0.05):
    assert (n-1) * sig_p > 1.0, 'p-value too small, increase number of samples'

    null_samples = []
    for _ in range(n):
        null_samples.append(np.sum(np.random.choice(scores, size=size)))
    return np.sort(null_samples, kind='mergesort')[int((1-sig_p)*n)]

# %%
@gin.configurable(denylist=['scores'])
def search(scores, sig_p=0.01, seed_size=2, max_size=20, extend_flanks=0):
    scores_cpy = copy.deepcopy(scores) # create copy for masking

    # dict of 
    kmer_sig_p_thresholds = {}

    discovered_kmers = []

    # from most to least important 2-mer
    for i in np.argsort(running_mean(scores, k=2)):
        # print('i:', i)
        extend_size = 0
        current_kmer_size = seed_size + 2*extend_size
        current_kmer_score = np.sum(scores_cpy[(i+extend_size):(i+2+extend_size)])
        current_kmer_sig = False
        while current_kmer_size <= max_size:
            # print(f'{current_kmer_score=}')
            # lazy compute significance thresholds for kmer size
            if current_kmer_size not in kmer_sig_p_thresholds:
                kmer_sig_p_thresholds[current_kmer_size] = find_significance_threshold(scores, size=current_kmer_size, sig_p=sig_p)

            # p = compute_p(scores, np.sum(scores[(i+extend_size):(i+2+extend_size)]), size=(size + 2*extend_size))
            # print(f'p={p}, range={(i-extend_size, i+extend_size+2)}')
            # check if kmer still significant
            if kmer_sig_p_thresholds[current_kmer_size] < current_kmer_score:
                current_kmer_sig = True
                extend_size += 1
                current_kmer_size = seed_size + 2*extend_size
                current_kmer_score = np.sum(scores_cpy[(i+extend_size):(i+2+extend_size)])
            else:
                break
        
        if current_kmer_sig: # only if kmer was significant
            kmer_start, kmer_stop = i-extend_size+1-extend_flanks, i+seed_size+extend_size+2-1+extend_flanks # get kmer range
            scores_cpy[kmer_start:kmer_stop] = -np.inf # mask kmer
            # print('size:', current_kmer_size)
            # print((kmer_start, kmer_stop))
            discovered_kmers.append((kmer_start, kmer_stop))

    return discovered_kmers