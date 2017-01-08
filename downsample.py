import pandas as pd
import os
import re

RAW_PATH = './reference_data'

lipid_match = re.compile("""
    (?P<hydroxy>[\d,]*OH-)?
    (?P<keto>\d*Oxo-)?
    (?P<cyclic>[\d,]*Cyc-)?
    (?P<methyl>[\d,]*Me-|Br-|A-|I-|Neo-)?
    (?P<base>C\d\d?:\d)
    (?P<db>ω?[\d,]+[ct]?)?
    (?P<other>-ALDE|-DMA|-DCA)?
""", re.VERBOSE)


def _colname_into_groups(colname):
    match = lipid_match.match(colname)
    if colname.startswith(('Other', 'Unknown')):
        return {'base': 'Other'}
    elif match is None:
        return None

    groups = {k: v for k, v in match.groupdict().items() if v is not None}
    # a little hacky, but force all ALDE/DMA/etc compounds into other
    if 'other' in groups:
        return {'base': 'Other'}
    return groups


def low_qual_name(colname):
    groups = _colname_into_groups(colname)
    if groups is None:
        return None
    name = ''

    if 'other' in groups:
        return 'Other'

    if 'hydroxy' in groups:
        name += 'OH-'
    if 'keto' in groups:
        name += 'Oxo-'
    if 'methyl' in groups:
        name += 'Br-'
    if 'cyclic' in groups:
        name += 'Cyc-'
    name += groups['base']
    name += groups.get('other', '')
    return name


def med_qual_name(colname):
    def get_pos(desc):
        return re.findall('\d+', desc)

    groups = _colname_into_groups(colname)
    if groups is None:
        return None
    name = ''

    for moiety in ['hydroxy', 'keto', 'cyclic']:
        if moiety in groups:
            if not get_pos(groups[moiety]):
                return 'Other'
            else:
                name += groups[moiety]

    if 'methyl' in groups:
        if groups['methyl'] == 'Br-':
            return 'Other'
        if groups['methyl'].endswith('Me-') and not get_pos(groups['methyl']):
            return 'Other'
        name += groups['methyl']

    name += groups['base']

    if not groups['base'].endswith(':0'):
        if 'db' not in groups:
            return 'Other'

        # only one double bond position is specified for medium quality
        # with more references, we could maybe bump this up to require
        # all the double bond positions (so medium quality is just
        # missing stereo information)
        name += 'ω' + str(get_pos(groups['db'])[0])

    return name


# TODO: quality is high if all of the stereochem is specified
# e.g. no db, or t/c specified (think this covers most cases?)
# need to check db/other group interactions though
# (also cyclic group sterochem is never annotated; assumed cis?)


lo_table = pd.DataFrame()
md_table = pd.DataFrame()
for filename in os.listdir(RAW_PATH):
    file_obj = open(os.path.join(RAW_PATH, filename), 'rb')
    raw_table = pd.read_csv(file_obj, sep='\t')

    # do a little cleanup; first, setting anything not filled into zero
    raw_table = raw_table.fillna(0)
    # now, set anything below the LOQ to be 1/3 the LOQ
    raw_table = raw_table.applymap(lambda x: float(x[1:]) / 3 if
                                   str(x).startswith('<') else x)

    # copy the data over into scratch tables to
    #   1. change column names by desired quality level
    #   2. merge now identical columns (e.g. Other's)
    for qual in ['lo', 'md']:
        table = pd.DataFrame()
        name_f = {
            'lo': low_qual_name,
            'md': med_qual_name,
        }[qual]
        for col in raw_table.columns:
            colname = name_f(col)
            if colname is None:
                # append * to the beginning of metadata columns
                table['*' + col] = raw_table[col]
            elif colname in table:
                table[colname] += raw_table[col].apply(float)
            else:
                table[colname] = raw_table[col].apply(float)

        if qual == 'lo':
            lo_table = pd.concat([lo_table, table])
        elif qual == 'md':
            md_table = pd.concat([md_table, table])

# replace all the NaN's with 0 or ""'s again (b/c of the
# concatenation
fillv = {c: ('' if c.startswith('*') else 0) for c in lo_table}
lo_table.fillna(fillv, inplace=True)
fillv = {c: ('' if c.startswith('*') else 0) for c in md_table}
md_table.fillna(fillv, inplace=True)

# now, quality filter out anything with more than 10% "Other"
lo_table = lo_table[lo_table['Other'] < 10]
md_table = md_table[md_table['Other'] < 10]

lo_table.to_csv('./lo_qual.tsv', sep='\t', index=False)
md_table.to_csv('./md_qual.tsv', sep='\t', index=False)
