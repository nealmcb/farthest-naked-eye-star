#!/usr/bin/env python3
"""
This script queries Gaia DR3 for stars with G < 6 and positive parallax,
computes two distance estimates (nominal and conservative lower-bound) in light years,
and then sorts the results to produce:
  - top 20 by nominal distance (expected distance)
  - top 20 by lower-bound distance
For each candidate the script uses SIMBAD to obtain extra fields:
  - SIMBAD main identifier (main_id)
  - A common name (if available)
  - Spectral type, from which a luminosity class is heuristically extracted
  - Variable star type (if present)
  - (Brightness range is not provided, so "N/A")
The final output tables use shortened headers:
  - main_id: SIMBAD main identifier
  - common: common name (if available)
  - src_id: Gaia DR3 source id
  - Gmag: Gaia phot_g_mean_mag (3 decimals)
  - plx: parallax (mas, 5 decimals)
  - plx_err: parallax_error (mas, 5 decimals)
  - d_nom: nominal distance in ly (rounded to integer)
  - d_lb: lower-bound distance in ly (rounded to integer)
  - spec: spectral type
  - lum: luminosity class
  - var: variable star type
  - brange: brightness range (always "N/A")
Two CSV files are written:
  - gaia_top20_nominal.csv for the nominal distance ranking
  - gaia_top20_lowerbound.csv for the conservative lower-bound ranking
"""

from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import re

###############################
# Part 1. Gaia Query and Distance Calculation
###############################

# Note: By default Gaia.launch_job returns at most 2000 rows.
# To get more rows, include a TOP clause in the query.
gaia_query = """
SELECT TOP 10000 source_id, ra, dec, phot_g_mean_mag, parallax, parallax_error
FROM gaiadr3.gaia_source
WHERE phot_g_mean_mag < 6
  AND parallax > 0
"""
print("Launching Gaia DR3 query …")
job = Gaia.launch_job(gaia_query)
gaia_results = job.get_results()
print("Gaia query returned {} rows.".format(len(gaia_results)))

# Conversion factor: 1 parsec = 3.26156 light years.
PC_TO_LY = 3.26156

def compute_distance_ly(parallax_mas):
    return (1000.0 / parallax_mas) * PC_TO_LY

def compute_lower_bound_distance_ly(parallax_mas, parallax_err_mas):
    return (1000.0 / (parallax_mas + parallax_err_mas)) * PC_TO_LY

nominal_distances = []
lower_bound_distances = []
for row in gaia_results:
    p = row['parallax']
    sigma = row['parallax_error']
    nominal_distances.append(compute_distance_ly(p))
    lower_bound_distances.append(compute_lower_bound_distance_ly(p, sigma))

gaia_results['distance_ly_nominal'] = nominal_distances
gaia_results['distance_ly_lower_bound'] = lower_bound_distances

# Create two sorted copies of the Gaia results:
sorted_nom = gaia_results.copy()
sorted_nom.sort('distance_ly_nominal', reverse=True)
top20_nom = sorted_nom[:20]

sorted_lb = gaia_results.copy()
sorted_lb.sort('distance_ly_lower_bound', reverse=True)
top20_lb = sorted_lb[:20]

###############################
# Part 2. SIMBAD Query Function and Helper Functions
###############################

custom_simbad = Simbad()
custom_simbad.reset_votable_fields()
# Request fields: MAIN_ID, ids, sp_type, and otype.
custom_simbad.add_votable_fields('MAIN_ID','ids','sp_type','otype')

def extract_luminosity_class(sp_type):
    """Extract luminosity class from the spectral type using a regex."""
    if sp_type is None or sp_type.strip() == "":
        return "N/A"
    m = re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
    return m.group(0) if m else "N/A"

def extract_common_name(ids_field, main_id):
    """
    Given the SIMBAD ids field (a string of alternative names separated by |) and the main_id,
    return a candidate common name that is not a Gaia identifier.
    """
    if ids_field is None:
        return "N/A"
    parts = [part.strip() for part in ids_field.split('|')]
    candidates = [x for x in parts if x != main_id]
    for candidate in candidates:
        if "gaia" in candidate.lower():
            continue
        # Heuristic: if candidate contains Greek letters or doesn't match a typical catalog pattern.
        if any(greek in candidate for greek in ['α','β','γ','δ','ε','ζ','η','θ','ι','κ','λ','μ','ν','ξ','ο','π','ρ','σ','τ','υ','φ','χ','ψ','ω']):
            return candidate
        if not re.match(r'^(HD|TYC|2MASS|HIP)\s*\d+', candidate, re.IGNORECASE):
            return candidate
    return "N/A"

def get_simbad_info(ra, dec, radius=2*u.arcsec):
    """
    Query SIMBAD by coordinate (with an increasing radius if needed) and return a dictionary with:
      - main_id
      - common_name
      - sp_type
      - lum_class
      - var_type
      - brightness_range (always "N/A")
    """
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    print(f"Querying SIMBAD at RA={ra:.5f}, Dec={dec:.5f}, radius={radius.to(u.arcsec)}")
    try:
        result = custom_simbad.query_region(coord, radius=radius)
    except Exception as e:
        print(f"SIMBAD query error at RA={ra:.5f}, Dec={dec:.5f}: {e}")
        result = None

    # Normalize column names to lowercase via rename_column.
    if result is not None:
        orig_cols = result.colnames.copy()
        for col in orig_cols:
            result.rename_column(col, col.lower())

    if result is None or len(result) == 0:
        if radius < 10*u.arcsec:
            return get_simbad_info(ra, dec, radius=radius + 3*u.arcsec)
        else:
            print(f"No SIMBAD match found for RA={ra:.5f}, Dec={dec:.5f} (up to radius {radius.to(u.arcsec)})")
            return {"main_id": "N/A", "common_name": "N/A", "sp_type": "N/A",
                    "lum_class": "N/A", "var_type": "N/A", "brightness_range": "N/A"}

    # Print diagnostic information for the first row.
    print("SIMBAD result row:")
    for col in result.colnames:
        print(f"  {col}: {result[col][0]}")
    
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
    common_name = extract_common_name(ids_field, main_id)
    
    return {"main_id": main_id,
            "common_name": common_name,
            "sp_type": sp_type if sp_type is not None else "N/A",
            "lum_class": lum_class,
            "var_type": var_type,
            "brightness_range": brightness_range}

###############################
# Part 3. Query SIMBAD for Each Top 20 Set
###############################

def enrich_with_simbad(top_table):
    sim_main_ids = []
    sim_common_names = []
    sim_sp_types = []
    sim_lum_classes = []
    sim_var_types = []
    sim_brightness_ranges = []
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
    # Add new columns
    top_table['simbad_main_id'] = sim_main_ids
    top_table['common_name'] = sim_common_names
    top_table['sp_type'] = sim_sp_types
    top_table['lum_class'] = sim_lum_classes
    top_table['var_type'] = sim_var_types
    top_table['brightness_range'] = sim_brightness_ranges
    return top_table

print("\nEnriching top20 nominal candidates with SIMBAD info …")
top20_nom = enrich_with_simbad(top20_nom)
print("\nEnriching top20 lower-bound candidates with SIMBAD info …")
top20_lb = enrich_with_simbad(top20_lb)

###############################
# Part 4. Build Final Tables with Shortened Headers and Format Numbers
###############################

def build_final_table(top_table):
    final_data = []
    for row in top_table:
        d_nom_int = int(round(row['distance_ly_nominal']))
        d_lb_int = int(round(row['distance_ly_lower_bound']))
        final_data.append({
            'main_id': row['simbad_main_id'],          # SIMBAD main identifier first
            'common': row['common_name'],              # then common name
            'src_id': row['source_id'],
            'Gmag': f"{row['phot_g_mean_mag']:.3f}",
            'plx': f"{row['parallax']:.5f}",
            'plx_err': f"{row['parallax_error']:.5f}",
            'd_nom': f"{d_nom_int}",
            'd_lb': f"{d_lb_int}",
            'spec': row['sp_type'],
            'lum': row['lum_class'],
            'var': row['var_type'],
            'brange': row['brightness_range']
        })
    # Create table with shortened headers.
    return Table(rows=final_data, names=['main_id','common','src_id','Gmag','plx','plx_err','d_nom','d_lb','spec','lum','var','brange'])

final_table_nom = build_final_table(top20_nom)
final_table_lb = build_final_table(top20_lb)

print("\nFinal Table (Top 20 Nominal Distance):")
print(final_table_nom)
print("\nFinal Table (Top 20 Lower-bound Distance):")
print(final_table_lb)

# Save the final tables to CSV files.
final_table_nom.write("gaia_top20_nominal.csv", format="csv", overwrite=True)
final_table_lb.write("gaia_top20_lowerbound.csv", format="csv", overwrite=True)
print("\nFinal tables saved to 'gaia_top20_nominal.csv' and 'gaia_top20_lowerbound.csv'.")
