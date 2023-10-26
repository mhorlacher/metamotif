# %%
import click
import gin
import tqdm
import numpy as np
from Bio import SeqIO

from metamotif.search import search

# %%
@click.command()
@click.argument('fasta', metavar='<sequences.fasta>')
@click.argument('scores', metavar='<scores.npy>') # help='A numpy array of shape (n_sequences, sequence_length).'
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