# %%
from pathlib import Path
from numpy import load

# %%
# example_data/QKI.attributions.ipy
QKI = load(str(Path(__file__).parent / 'QKI.attributions.npy'))