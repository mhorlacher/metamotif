# %%
import argparse
from pathlib import Path

from tqdm import tqdm

from metamotif.alignment import VariableLengthSeededMotifAlignment
from metamotif.utils import sequence2onehot, write_motif_tsv
from metamotif.visualize import plot_motif

# %%
def load_kmers(tsv, to_onehot=True):
    kmers = []
    with open(tsv) as f:
        for i, line in enumerate(f):
            name, kmer, score, *_ = line.strip().split('\t')
            if len(kmer) < 2:
                print('-->', i, kmer)
                exit()
            if to_onehot:
                kmer = sequence2onehot(kmer)
            kmers.append((kmer, score))
    return kmers

def find_motifs(kmers, min_agreement = 3, min_frac_agreement = .5):
    seeded_alignments = [VariableLengthSeededMotifAlignment(kmers[0])]
    for kmer in tqdm(kmers[1:], total=len(kmers)-1):
        if len(kmer) < 2:
            print(kmer)
        for alignment in seeded_alignments:
            score_threshold = max(min_agreement, min_frac_agreement * min(len(kmer), len(alignment.seed)))
            if alignment.align(kmer, min_score = score_threshold):
                break
        else:
            seeded_alignments.append(VariableLengthSeededMotifAlignment(kmer))
    return seeded_alignments

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('kmer_csv')
    parser.add_argument('--min-support', type=int, default=100)
    parser.add_argument('--max-motifs', type=int, default=5)
    parser.add_argument('--min-alignment-agreement', type=int, default=3)
    parser.add_argument('--min-alignment-agreement-frac', type=float, default=.5)
    parser.add_argument('-o', '--output-directory')
    args = parser.parse_args()

    # set total support
    total_support = len(load_kmers(args.kmer_csv, to_onehot=False))
    
    # load kmers
    kmers = load_kmers(args.kmer_csv)
    kmers = sorted(kmers, key = lambda x: x[1], reverse=True)
    kmers = [kmer for kmer, score in kmers]
    
    # find motifs
    motifs = find_motifs(kmers, min_agreement = args.min_alignment_agreement, min_frac_agreement = args.min_alignment_agreement_frac)
    
    # save/plot motifs
    output_path = Path(args.output_directory) / 'motif-{i}'
    output_path.parent.mkdir(exist_ok=True)
    output_path = str(output_path)
    for i, motif in enumerate(sorted(motifs, key = lambda x: x.support, reverse=True)):
        if i >= args.max_motifs:
            break
        if motif.support < args.min_support:
            break
        
        write_motif_tsv(motif.pwm, filepath=(output_path.format(i=i) + '.tsv'), meta_info={'support': motif.support})
        fig = plot_motif(motif.pwm, ylab = 'Occupancy', title=f'{motif.support}/{total_support}') # {len(kmers)}
        fig.savefig((output_path.format(i=i) + '.pdf'), bbox_inches='tight')
        fig.savefig((output_path.format(i=i) + '.png'), bbox_inches='tight')

# %%
if __name__ == '__main__':
    main()