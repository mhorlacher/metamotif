# %%
import logging
import copy

import numpy as np

from metamotif.utils import pad_matrix

# %%
class SeededMotifAlignment:
    # def __init__(self, seed, size=5, min_agreement=3):
    def __init__(self, seed, min_agreement=3):
        # params
        self.size = len(seed)
        self.min_agreement = min_agreement
        self.padding = self.size - self.min_agreement
        self.extended_size = self.size + 2 * self.padding
        
        # seed
        self.seed = seed
        self.padded_seed = np.zeros(shape=(self.extended_size, 4), dtype=np.float32)
        self.padded_seed[self.padding:(-self.padding)] = self.seed
        
        # pam
        self.pam = np.zeros(shape=(self.extended_size, 4), dtype=np.float32)
        self.pam[self.padding:(-self.padding)] = self.seed
        # self.support = np.zeros(shape=(self.extended_size, ))
        # self.support[self.padding:(-self.padding)] = np.ones(shape=(self.size, ))
        
        self.support = 1
    
    @property
    def pwm(self):
        return self.pam / self.support
    
    def align(self, kmer):
        assert len(kmer) == self.size
        
        agreement_max_idx, agreement_max_val = None, 0
        for i in range(self.extended_size - self.size + 1):
            agreement_i_val = np.sum(kmer * self.padded_seed[i:(i+self.size)])
            if agreement_i_val > agreement_max_val:
                agreement_max_val = agreement_i_val
                agreement_max_idx = i
                
        if agreement_max_val < self.min_agreement:
            return False
        else:
            self.pam[agreement_max_idx:(agreement_max_idx + self.size)] += kmer
            self.support += 1
            return True

# %%
class VariableLengthSeededMotifAlignment:
    def __init__(self, seed):
        # seed
        self.seed = copy.deepcopy(seed)
        
        # alignment
        self.pfm = self.seed
        self.pfm_seed_offset = 0
        self.support = 1
    
    @property
    def pwm(self):
        return self.pfm / self.support
    
    def align(self, motif, sim_fn = lambda x, y: np.sum(x * y), min_score = None):
        motif = copy.deepcopy(motif)
        assert len(motif) > 1

        # assign short and long motif
        motif_long, motif_short = (motif, self.seed) if len(motif) > len(self.seed) else (self.seed, motif)
        
        # add padding to long motif
        padding_length =  len(motif_short) - 1 # hard-coded min overlap of 1
        motif_long = pad_matrix(motif_long, padding_left=padding_length, padding_right=padding_length)

        # scan motif_long for best alignment of motif_short
        max_score = -np.inf
        max_score_offset = None
        for i in range(len(motif_long) - len(motif_short) + 1):
            # take subsequence from motif_long (of length of motif_short) at offset i
            motif_long_subseq = motif_long[i:(i+len(motif_short))]

            # compute alignment score
            alignment_score = sim_fn(motif_long_subseq, motif_short)

            if alignment_score > max_score:
                max_score = alignment_score
                max_score_offset = i
        
        # break if score is too low
        if min_score is not None:
            if max_score < min_score:
                return False
        
        # align motif_short to motif_long
        logging.debug(f'{max_score_offset=}, {max_score=}')
        motif_aligned = motif_long
        motif_aligned[max_score_offset:(max_score_offset+len(motif_short))] += motif_short

        # now, we can find the offset of the aligned motif relative to the seed
        if len(motif) > len(self.seed):
            # if the non-seed motif is the longer sequence
            offset_to_seed = -max_score_offset + padding_length
        else:
            # if the non-seed motif is the shorter sequence, it's the offset of the alignment, minus the padding
            offset_to_seed = max_score_offset - padding_length

        
        # --> using the offset, we can add the aligned motif to the PFM
        logging.debug(f'{self.pfm_seed_offset=}, {offset_to_seed=}, {len(motif)=}, {len(self.pfm)=}')

        # first, we need to check if we need to left-pad the current PFM (if the motif starts BEFORE the seed)
        if (self.pfm_seed_offset + offset_to_seed) < 0:
            logging.debug(f'pad-left: {-(self.pfm_seed_offset + offset_to_seed)}')
            self.pfm = pad_matrix(self.pfm, padding_left = -(self.pfm_seed_offset + offset_to_seed))
            # ... and update the offset!
            self.pfm_seed_offset += -(self.pfm_seed_offset + offset_to_seed)
        
        # next, we need to check if we need to right-pad the current PFM
        if (len(self.pfm) - self.pfm_seed_offset - offset_to_seed) < len(motif):
            logging.debug(f'pad-right: {len(motif) - (len(self.pfm) - self.pfm_seed_offset - offset_to_seed)}')
            # in this case, we need to right-pad the current PFM
            self.pfm = pad_matrix(self.pfm, padding_right = len(motif) - (len(self.pfm) - self.pfm_seed_offset - offset_to_seed))

        # finally, we can add the aligned motif to the PFM (using the offset of the seed and the offset of the motive relative to the seed)
        self.pfm[(self.pfm_seed_offset + offset_to_seed):(self.pfm_seed_offset + offset_to_seed + len(motif))] += motif
        self.support += 1

        return max_score