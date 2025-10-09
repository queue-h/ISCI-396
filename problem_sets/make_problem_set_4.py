import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table

# %%

table_1 = Table.read(
    '/mnt/backup/Catalogs/COSMOS/2020/COSMOS2020_FARMER_R1_v2.2_p3.fits')
table_1 = table_1[(table_1['FLAG_COMBINED'] == 0) &
                  np.isin(table_1['lp_type'], [0, 1, 2])]
table_2 = Table.read(
    '/mnt/backup/Catalogs/DESI/COSMOS-XMM/DESI-COSMOS-v2.0.fits')
table_2 = table_2[table_2['QUALITY_Z']]

# %%


def match(table_1, table_2, match_radius=1 * u.arcsec):
    """Match two catalogs spatially.

    Parameters
    ----------
    table_1 : astropy.table.Table
        Catalog that matches are searched for. Must have columns `RA` and
        `DEC`.
    table_2 : astropy.table.Table
        Catalog in which matches for the first table are looked for. Must have
        columns `RA` and `DEC`.
    match_radius : astropy.units.quantity.Quantity, optional
        Matching radius required for a succesful match. Default is 1.0 arcsec.

    Returns
    -------
    match : numpy.ndarray
        Array with the same length as `table_1` indicating whether a match
        within the matching radius was found.
    idx : numpy.ndarray
        Array with the same length as `table_1` storing the index of the
        nearest match in `table_2`.

    """
    c_1 = SkyCoord(table_1['RA'], table_1['DEC'], unit='deg')
    c_2 = SkyCoord(table_2['RA'], table_2['DEC'], unit='deg')
    idx, sep2d = match_coordinates_sky(c_1, c_2)[:2]
    return sep2d < match_radius, idx


table_1.rename_columns(['ALPHA_J2000', 'DELTA_J2000'], ['RA', 'DEC'])
table_2.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])

match, idx = match(table_1, table_2)
idx = idx[match]
table_1 = table_1[match]

# %%

features = []

for key in table_1.colnames:
    if key[-5:] == '_FLUX' and key not in [
            'GALEX_FUV_FLUX', 'UVISTA_NB118_FLUX']:
        features.append(key)

table_1['Z_SPEC'] = table_2['Z'][idx]
table_1['Z_PHOT'] = np.where(table_1['lp_type'] != 2, table_1['lp_zBEST'],
                             table_1['lp_zq'])
table_1.keep_columns(features + ['Z_PHOT', 'Z_SPEC'])

# %%

mask = np.zeros(len(table_1), dtype=bool)
for key in table_1.colnames:
    try:
        mask = mask | table_1[key].mask
    except AttributeError:
        pass
table_1 = table_1[~mask]
table_1.write('Problem_Set_4_Redshifts.csv', overwrite=True)
