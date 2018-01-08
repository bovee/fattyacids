from math import log

import pandas as pd
from ete3 import NodeStyle, TextFace, Tree, TreeStyle, NCBITaxa

ncbi = NCBITaxa()

l = pd.read_csv('lo_qual.tsv', sep='\t', quoting=3)
new_l = l[['*Tax ID']]
new_l = new_l.assign(br=0, oh=0, cyc=0, long=0, unsat=0)

def get_genus(tax_id):
    genus = [i for i in ncbi.get_lineage(tax_id) if ncbi.get_rank([i])[i] == 'genus']
    if len(genus) > 0:
        return str(genus[0])
    else:
        return None

def logt(x):
    return ((log(x + 0.01) / log(100)) + 1) * 50

for col in l.columns:
    if col.startswith('Br-') or '-Br-' in col:
        new_l['br'] += l[col]
    if col.startswith('OH-'):
        new_l['oh'] += l[col]
    if col.startswith('Cyc-') or '-Cyc-' in col:
        new_l['cyc'] += l[col]
    if 'C2' in col or 'C3' in col:
        new_l['long'] += l[col]
    if col.endswith(':0'):
        new_l['unsat'] += l[col]

new_l['genus'] = [get_genus(i) for i in new_l['*Tax ID']]
lookup = {r['genus']: dict(r) for _, r in new_l.groupby('genus').mean().reset_index().iterrows()}

hug_tree = Tree(HUG_TREE)

shared_genus = set(hug_tree.get_leaf_names()).intersection(new_l['genus'])
hug_tree.prune(shared_genus)

ts = TreeStyle()
ts.show_leaf_name = False
ts.optimal_scale_level = 'full'
ts.min_leaf_separation = 2
ts.mode = "c"

for n in hug_tree.traverse():
    nstyle = NodeStyle()
    if n.is_leaf():
        r = logt(lookup[n.name]['br'])
        g = lookup[n.name]['unsat']
        b = lookup[n.name]['long']
        bright_color = '#{0:02x}{1:02x}{2:02x}'.format(
            int(r * 2.55),
            int(g * 2.55),
            int(b * 2.55)
        )
        color = '#{0:02x}{1:02x}{2:02x}'.format(
            int(r) + 155,
            int(g) + 155,
            int(b) + 155
        )
        nstyle['fgcolor'] = bright_color
        nstyle['bgcolor'] = color
        nstyle['hz_line_color'] = '#404040'
        nstyle['size'] = 200

        getattr(n.faces, 'branch-right')[0] = []
        t = TextFace(
            text=ncbi.get_taxid_translator([n.name])[int(n.name)],
            fsize=240,
            ftype='Verdana',
            penwidth=1,
        )
        n.add_face(t, 0)
    else:
        nstyle['vt_line_color'] = '#404040'
        nstyle['hz_line_color'] = '#404040'
        nstyle['size'] = 0
    n.set_style(nstyle)

hug_tree.show(tree_style=ts)
