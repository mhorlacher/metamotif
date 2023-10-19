# %%
#import argh
import click

# %%
from metamotif.bin import search

# %%
@click.group()
def main():
    pass

# %%
main.add_command(search.main, name='search')
# %%
if __name__ == '__main__':
    main()