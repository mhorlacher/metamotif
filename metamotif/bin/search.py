# %%
import copy

import click
import gin
import tqdm
import numpy as np
from Bio import SeqIO

from metamotif.utils import onehot2sequence, running_mean

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

# %%
@click.command()
@click.argument('fasta', metavar='<sequences.fasta>')
@click.argument('scores', metavar='<scores.npy>')
@click.option('-c', '--config', default=None)
@click.option('-o', '--output', default=None)
# @click.option('--alphabet', default='ACGT')
def main(fasta, scores, config, output):
    if config is not None:
        gin.parse_config_file(config)
    
    with open(output, 'w') as fout:
        # for i, (sequence_oh, scores_oh) in enumerate(np.load(data)):
        for record, scores_ in tqdm.tqdm(zip(SeqIO.parse(fasta, 'fasta'), np.load(scores))):
            # print(scores_.shape)
            for (kmer_start, kmer_stop) in search(scores_):
                kmer_seq = str(record.seq[kmer_start:kmer_stop])
                kmer_score = np.sum(scores_[kmer_start:kmer_stop])
                
                print(f'{record.id}\t{kmer_seq}\t{kmer_score:.4f}\t{kmer_start}\t{kmer_stop}\t{kmer_stop-kmer_start}', file=fout)

# %%
if __name__ == '__main__':
    main()