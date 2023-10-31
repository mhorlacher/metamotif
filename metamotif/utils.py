# %%
import numpy as np

# %%
base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def sequence2int(sequence):
    return [base2int.get(base, 999) for base in sequence]

def sequence2onehot(sequence, sigma='ACGT'):
    try:
        import tensorflow as tf
    except ImportError:
        raise ImportError('tensorflow is required for onehot encoding')

    _sigma_to_int = dict(zip(sigma, range(len(sigma))))
    return tf.one_hot([_sigma_to_int.get(x, -1) for x in sequence], depth=len(sigma)).numpy()

def onehot2sequence(onehot, sigma='ACGT'):
    assert onehot.shape[1] == len(sigma)
    return ''.join([(sigma[np.argmax(col)]) if np.max(col) > 0 else '0' for col in onehot])

# %%
def running_mean(x, k):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[k:] - cumsum[:-k]) / float(k)

# %%
def matrix_to_transfac(mtrx, id=None, alphabet='ACGT'):
    assert len(alphabet) == mtrx.shape[1]
    mtrx = np.array(mtrx, dtype=np.int64)

    transfac = ''
    transfac += f'AC {id}\n' if id is not None else 'ID\n'
    transfac += f'ID {id}\n' if id is not None else 'ID\n'
    transfac += 'P0\t' + '\t'.join(alphabet) + '\n'
    for i in range(mtrx.shape[0]):
        transfac += str(i+1).zfill(2) + '\t' + '\t'.join([str(x) for x in mtrx[i]]) + '\n'
    transfac += 'XX\n'
    transfac += '//\n'

    return transfac

# %%
def matrix_to_meme(mtrx, id=None, alphabet='ACGT'):
    raise NotImplementedError

# %%
def transfac_to_matrix(transfac_string):
    raise NotImplementedError

# %%
def write_motif_tsv(motif_array, filepath, sigma=['A', 'C', 'G', 'U'], meta_info={}):
    assert motif_array.shape[1] == len(sigma)
    with open(filepath, 'w') as f:
        for key, value in meta_info.items():
            print(f'#{key}={value}', file=f)
        print('\t'.join(sigma), file=f)
        for row in motif_array:
            print('\t'.join(map(str, row)), file=f)

# %%
def pad_matrix(matrix, padding_left=0, padding_right=0):
    """Pads the input position-weight matrix. 

    Args:
        matrix (tf.Tensor): Input matrix.
        padding (int, optional): Padding size. Defaults to 0.

    Returns:
        tf.Tensor: Padded PWM. 
    """

    try:
        import tensorflow as tf
    except ImportError:
        raise ImportError('tensorflow is required for padding')

    return tf.pad(matrix, [[padding_left, padding_right,], [0, 0]], 'CONSTANT').numpy()

def remove_padding(matrix, padding_value=0):
    """Removes padding from the input matrix. 

    Args:
        matrix (tf.Tensor): Input matrix.
        padding_value (int, optional): Padding value. Defaults to 0.

    Returns:
        tf.Tensor: Padded matrix. 
    """
    matrix = np.array(matrix)
    return matrix[np.sum(matrix, axis=1) != padding_value]