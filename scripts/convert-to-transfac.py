import numpy as np
import click

from metamotif.utils import sequence2onehot, matrix_to_transfac

# %%
@click.command()
@click.argument('kmers', metavar='<kmers.tsv>', required=True)
@click.option('--min-count', type=int, default=2)
@click.option('--alphabet', default='ACGT')
def main(kmers, min_count, alphabet):
    with open(kmers) as f:
        kmer_counts = {}
        for line in f:
            kmer = line.strip().split('\t')[1]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1
    
    i = 0
    for kmer, count in sorted(kmer_counts.items(), key = lambda x: x[1], reverse=True):
        if count < min_count:
            continue
        kmer_onehot = sequence2onehot(kmer, sigma=alphabet)
        kmer_onehot *= count
        print(matrix_to_transfac(kmer_onehot, id=i, alphabet=alphabet), end='')
        i += 1

# %%
if __name__ == '__main__':
    main()