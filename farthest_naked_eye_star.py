#!/usr/bin/env python3
"""
Find the most distant stars visible to the naked eye.

Usage
   python farthest_naked_eye_star.py [query]

Produces two csv files based on nominal parallax, and a lower bound
by adding the parallax error bound to the nominal value:

* gaia_top20_nominal.csv
* gaia_top20_lowerbound.csv

If a query is provided, all the matching stars will be described.
For example,
`python farthest_naked_eye_star.py "source_id in (1998148532777850880, 5350358584482202880)"`
returns a table of data for just two stars, Rho Cas and Eta Carinae, ignoring the fact
that the latter doesn't have a parallax in Gaia DR3.

This script implements the following strategy:
  - Query Gaia DR3 for stars with Gmag < 6.2 and positive parallax.
  - Compute nominal and conservative lower-bound distances (in ly).
  - Form two candidate tables: one sorted by nominal distance and one by lower-bound distance.
  - Initially select the top 50 candidates from each.
  - Enrich each candidate with SIMBAD data:
      * SIMBAD main identifier, common name, Vmag, spectral type, luminosity class, variable star type.
  - Then filter each table to keep only those with Vmag < 6.
  - Finally, build and save final tables with the following columns:
      main_id, common, Gaia DR3, Vmag, Gmag, ly_nom, ly_lb, spec, lum, var, brange, Con, RA, Dec.
  - Print how many objects remain in each selection.
"""

import sys
import logging
import numpy as np
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import re

# logging.basicConfig(level=logging.INFO)

vmag_limit = 6.0
show_all = False
where = None

###############################
# Part 1. Gaia Query and Distance Calculation
###############################

# Query Gaia DR3 for stars with Gmag < 6.2, ignoring those with invalid parallax
# Note, this Gaia DR3 query returned 8318 rows.
gaia_query = """
SELECT TOP 10000 source_id, ra, dec, phot_g_mean_mag, parallax, parallax_error
FROM gaiadr3.gaia_source
WHERE phot_g_mean_mag < 6.2
"""

# Given a WHERE clause, query for just that
if len(sys.argv) > 1:
    where = sys.argv[1]
    print(f'Custom query {where=}')
    gaia_query = f"""
      SELECT source_id, ra, dec, phot_g_mean_mag, parallax, parallax_error
      FROM gaiadr3.gaia_source
      WHERE {where}"""
    show_all = True


print("Launching Gaia DR3 query …")
job = Gaia.launch_job(gaia_query)
gaia_results = job.get_results()
print(f"Gaia query returned {len(gaia_results)} rows.")

# Conversion factor: 1 parsec = 3.26156 light years.
PC_TO_LY = 3.26156

def compute_distance_ly(parallax_mas):
    "Given parallax in mas, compute distance in light years"

    return (1000.0 / parallax_mas) * PC_TO_LY

def compute_lower_bound_distance_ly(parallax_mas, parallax_err_mas):
    "Given parallax in mas, increased by STD, return reasonable lower-bound distance in light years"

    return (1000.0 / (parallax_mas + parallax_err_mas)) * PC_TO_LY

nominal_distances = []
lower_bound_distances = []
for row in gaia_results:
    plx = row['parallax']
    sigma = row['parallax_error']
    if plx < 0:
        logging.info(f"parallax < 0: {plx=}, {sigma=}, {row}")
        d = float('NaN')
    else:
        d = compute_distance_ly(plx)
    nominal_distances.append(d)
    if plx + sigma < 0:
        logging.info(f"lower bound < 0: plx + sigma < 0 {plx=}, {sigma=}, {row}")
        d = float('NaN')
    else:
        d = compute_lower_bound_distance_ly(plx, sigma)
    lower_bound_distances.append(d)

gaia_results['distance_ly_nominal'] = nominal_distances
gaia_results['distance_ly_lower_bound'] = lower_bound_distances

# Unless we're reporting on everything found for a custom query, make two sorted copies:
# FIXME: For custom queries, either sort while including unsortable rows also,
#  or leave unsorted but only process one table
if not where:
    sorted_nom = gaia_results[gaia_results['distance_ly_nominal'] > 0].copy()
    sorted_nom.sort('distance_ly_nominal', reverse=True)
    top50_nom = sorted_nom[:50]

    sorted_lb = gaia_results[gaia_results['distance_ly_lower_bound'] > 0].copy()
    sorted_lb.sort('distance_ly_lower_bound', reverse=True)
    top50_lb = sorted_lb[:50]
else:
    top50_nom = gaia_results
    top50_lb = gaia_results

###############################
# Part 2. SIMBAD Query and Helper Functions
###############################

custom_simbad = Simbad()
custom_simbad.reset_votable_fields()
# Request: MAIN_ID, ids, sp_type, otype, and visual magnitude ("v").
custom_simbad.add_votable_fields('MAIN_ID','ids','sp_type','otype','V')

def extract_luminosity_class(sp_type):
    if sp_type is None or sp_type.strip() == "":
        return "N/A"
    m = re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
    return m.group(0) if m else "N/A"

def extract_common_name(ids_field, main_id):
    if ids_field is None:
        return "N/A"
    parts = [part.strip() for part in ids_field.split('|')]
    candidates = [x for x in parts if x != main_id and "gaia" not in x.lower()]
    for candidate in candidates:
        if any(greek in candidate for greek in ['α','β','γ','δ','ε','ζ','η','θ','ι','κ','λ','μ','ν','ξ','ο','π','ρ','σ','τ','υ','φ','χ','ψ','ω']):
            return candidate
        if not re.match(r'^(HD|TYC|2MASS|HIP)\s*\d+', candidate, re.IGNORECASE):
            return candidate
    return "N/A"

def get_simbad_info(ra, dec, radius=2*u.arcsec):
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    logging.info(f"Querying SIMBAD at RA={ra:.5f}, Dec={dec:.5f}, radius={radius.to(u.arcsec)}")
    try:
        result = custom_simbad.query_region(coord, radius=radius)
    except Exception as e:
        logging.error(f"SIMBAD query error at RA={ra:.5f}, Dec={dec:.5f}: {e}")
        result = None

    # Normalize column names to lowercase.
    if len(result) > 0:
        orig_cols = result.colnames.copy()
        for col in orig_cols:
            result.rename_column(col, col.lower())

    if result is None or len(result) == 0:
        if radius < 10*u.arcsec:
            return get_simbad_info(ra, dec, radius=radius + 3*u.arcsec)
        else:
            logging.error(f"No SIMBAD match found for RA={ra:.5f}, Dec={dec:.5f} (up to radius {radius.to(u.arcsec)})")
            return {"main_id": "N/A", "common_name": "N/A", "sp_type": "N/A",
                    "lum_class": "N/A", "var_type": "N/A", "brightness_range": "N/A", "vmag": "N/A"}

    logging.info("SIMBAD result row:")
    for col in result.colnames:
        logging.info(f"  {col}: {result[col][0]}")
    
    main_id = result['main_id'][0] if 'main_id' in result.colnames else "N/A"
    if isinstance(main_id, bytes):
        main_id = main_id.decode('utf-8')
        
    ids_field = result['ids'][0] if 'ids' in result.colnames else "N/A"
    if isinstance(ids_field, bytes):
        ids_field = ids_field.decode('utf-8')
        
    sp_type = result['sp_type'][0] if 'sp_type' in result.colnames and result['sp_type'][0] is not None else "N/A"
    if isinstance(sp_type, bytes):
        sp_type = sp_type.decode('utf-8')
        
    lum_class = extract_luminosity_class(sp_type) if sp_type != "N/A" else "N/A"
    
    otype = result['otype'][0] if 'otype' in result.colnames and result['otype'][0] is not None else "N/A"
    if isinstance(otype, bytes):
        otype = otype.decode('utf-8')
    var_type = otype if otype is not None and ("Var" in otype or "V*" in otype) else "N/A"
    
    brightness_range = "N/A"  # Not provided by SIMBAD.

    if 'v' in result.colnames and result['v'][0] is not None:
        vmag = result['v'][0]
        if isinstance(vmag, bytes):
            vmag = vmag.decode('utf-8')
    else:
        vmag = "N/A"
    
    common_name = extract_common_name(ids_field, main_id)
    
    return {"main_id": main_id,
            "common_name": common_name,
            "sp_type": sp_type if sp_type is not None else "N/A",
            "lum_class": lum_class,
            "var_type": var_type,
            "brightness_range": brightness_range,
            "vmag": vmag}

###############################
# Part 3. Enrich the Candidate Tables with SIMBAD Data
###############################

def enrich_with_simbad(top_table):
    """Add SIMBAD data to top_table
    FIXME: this seems overly verbose
    """

    sim_main_ids = []
    sim_common_names = []
    sim_sp_types = []
    sim_lum_classes = []
    sim_var_types = []
    sim_brightness_ranges = []
    sim_vmags = []
    for row in top_table:
        ra = row['ra']
        dec = row['dec']
        sim_info = get_simbad_info(ra, dec)
        sim_main_ids.append(sim_info["main_id"])
        sim_common_names.append(sim_info["common_name"])
        sim_sp_types.append(sim_info["sp_type"])
        sim_lum_classes.append(sim_info["lum_class"])
        sim_var_types.append(sim_info["var_type"])
        sim_brightness_ranges.append(sim_info["brightness_range"])
        sim_vmags.append(sim_info["vmag"])
    top_table['simbad_main_id'] = sim_main_ids
    top_table['common_name'] = sim_common_names
    top_table['sp_type'] = sim_sp_types
    top_table['lum_class'] = sim_lum_classes
    top_table['var_type'] = sim_var_types
    top_table['brightness_range'] = sim_brightness_ranges
    top_table['vmag'] = sim_vmags
    return top_table

print("\nEnriching nominal candidates (top50_nom) with SIMBAD info …")
top50_nom = enrich_with_simbad(top50_nom)

print("\nEnriching lower-bound candidates (top50_lb) with SIMBAD info …")
top50_lb = enrich_with_simbad(top50_lb)

def filter_by_vmag(top_table, show_all):
    "Filter table to keep only rows with Vmag < vmag_limit unless show_all"

    filtered_rows = []
    for row in top_table:
        try:
            if isinstance(row['vmag'], np.ma.core.MaskedArray):
                vmag = float('NaN')
                logging.info(f"Bypassing {row['source_id']=} vmag MaskedArray")
            else:
                vmag = float(row['vmag'])
        except (ValueError, TypeError):
            logging.error(f'vmag ValueError in {row=}')
            continue
        if vmag < vmag_limit or show_all:
            filtered_rows.append(row)

    assert len(filtered_rows) > 0, f'No rows left after filtering: {top_table}'
    return Table(rows=filtered_rows, names=top_table.colnames)

# Now, filter each table to keep only rows with Vmag < vmag_limit.

filtered_nom = filter_by_vmag(top50_nom, show_all)
filtered_lb = filter_by_vmag(top50_lb, show_all)

print(f"\nNumber of nominal candidates: {len(filtered_nom)}")
print(f"Number of lower-bound candidates: {len(filtered_lb)}")

###############################
# Part 4. Build Final Tables with New Column Order and Labels
###############################

def build_final_table(top_table):
    final_data = []
    for row in top_table:
        try:
            d_lb_int = int(round(row['distance_ly_lower_bound']))
        except:
            logging.info(f"Error distance_ly_lower_bound non-int {row}")
            d_lb_int = -1
        try:
            d_nom_int = int(round(row['distance_ly_nominal']))
        except:
            logging.info(f"Error distance_ly_nominal non-int {row}")
            d_nom_int = -1

        # Compute constellation abbreviation
        coord = SkyCoord(ra=row['ra']*u.deg, dec=row['dec']*u.deg, frame='icrs')
        try:
            constellation = coord.get_constellation(short_name=True)
        except TypeError:
            constellation = coord.get_constellation()

        # Avoid pesky elusive FutureWarning: Format strings passed to MaskedConstant are ignored, but in future may error or produce different behavior
        # Note, masked value can be accessed as row['parallax'].data, typically 0.0 here
        #  Or coerced to NaN via float(row['parallax'])
        if isinstance(row['parallax'], np.ma.core.MaskedArray):
            logging.info(f"Fixing {row['source_id']=} parallax MaskedArray formatting")
            row['parallax'] = float('NaN')
            row['parallax_error'] = float('NaN')
        final_data.append({
            'main_id': row['simbad_main_id'],         # SIMBAD main identifier
            'common': row['common_name'],             # common name
            'Vmag': f"{float(row['vmag']):.3f}" if row['vmag'] != "N/A" else "N/A",
            'Gmag': f"{row['phot_g_mean_mag']:.3f}",
            'plx': f"{row['parallax']:.5f}",           # Parallax in mas
            'plx_err': f"{row['parallax_error']:.5f}",  # Parallax error in mas
            'ly_nom': f"{d_nom_int}",
            'ly_lb': f"{d_lb_int}",
            'Con': constellation,
            'RA': f"{row['ra']:.5f}",
            'Dec': f"{row['dec']:.5f}",
            'spec': row['sp_type'],
            'lum': row['lum_class'],
            'var': row['var_type'],
            'Gaia DR3': row['source_id']             # Gaia source id
        })
    return Table(rows=final_data)

final_table_nom = build_final_table(filtered_nom)
final_table_lb = build_final_table(filtered_lb)

print("\nFinal Table (Nominal candidates):")
print(final_table_nom)
print("\nFinal Table (Lower-bound candidates):")
print(final_table_lb)

# Save to CSV:
final_table_nom.write("gaia_top20_nominal.csv", format="csv", overwrite=True)
final_table_lb.write("gaia_top20_lowerbound.csv", format="csv", overwrite=True)

print("\nFinal tables saved to 'gaia_top20_nominal.csv' and 'gaia_top20_lowerbound.csv'.")
