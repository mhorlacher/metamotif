# %%
import argparse
from pathlib import Path

from tqdm import tqdm

from metamotif.alignments import SeededMotifAlignment
from metamotif.utils import sequence2onehot, write_motif_tsv
from metamotif.visualize import plot_motif

# %%
def load_kmers(tsv, to_onehot=True):
    kmers = []
    with open(tsv) as f:
        _ = f.readline()
        for line in f:
            name, kmer, score = line.strip().split('\t')
            if to_onehot:
                kmer = sequence2onehot(kmer)
            kmers.append((kmer, score))
    return kmers

def find_motifs(kmers):
    seeded_alignments = [SeededMotifAlignment(kmers[0][0])]
    for kmer, score in tqdm(kmers, total=len(kmers)):
        for alignment in seeded_alignments:
            success = alignment.align(kmer)
            if success:
                break
        else:
            seeded_alignments.append(SeededMotifAlignment(kmer))
    return seeded_alignments

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('kmer_csv')
    parser.add_argument('--min-support', type=int, default=100)
    parser.add_argument('--max-motifs', type=int, default=5)
    parser.add_argument('-o', '--output-directory')
    args = parser.parse_args()
    
    # load kmers
    kmers = load_kmers(args.kmer_csv)
    kmers = sorted(kmers, key = lambda x: x[1], reverse=True)
    
    # find motifs
    motifs = find_motifs(kmers)
    
    # save/plot motifs
    output = Path(args.output_directory)
    output.mkdir(exist_ok=True)
    for i, motif in enumerate(sorted(motifs, key = lambda x: x.support, reverse=True)):
        if i >= args.max_motifs:
            break
        if motif.support < args.min_support:
            break
        
        write_motif_tsv(motif.pwm, filepath=(output / f'motif-{i}.tsv'), meta_info={'support': motif.support})
        fig = plot_motif(motif.pwm, ylab = 'Occupancy', title=f'{motif.support}/5000') # {len(kmers)}
        fig.savefig(output / f'motif-{i}.pdf', bbox_inches='tight')

# %%
if __name__ == '__main__':
    main()