#!/usr/bin/env python3
"""
Query Gaia DR3 for stars with G < 6, compute distances (in light years)
from the parallax (and from parallax+error for a conservative lower bound),
then for the top 20 objects retrieve SIMBAD name and star–type info.
The final table columns (in order) are:
  - Common Name (if available)
  - SIMBAD Main Identifier
  - Gaia DR3 source_id
  - Gmag
  - parallax (mas)
  - parallax_error (mas)
  - nominal distance (ly)
  - lower-bound distance (ly)
  - Spectral Type
  - Luminosity Class
  - Variable Star Type
  - Brightness Range

Numeric columns are formatted with three decimals for magnitude,
and five for parallaxes and distances.
"""

from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import re

###############################
# Part 1. Gaia Query and Distance Calculation
###############################

query = """
SELECT source_id, ra, dec, phot_g_mean_mag, parallax, parallax_error
FROM gaiadr3.gaia_source
WHERE phot_g_mean_mag < 6
  AND parallax > 0
"""
print("Launching Gaia DR3 query …")
job = Gaia.launch_job(query)
gaia_results = job.get_results()
print("Gaia query returned {} rows.".format(len(gaia_results)))

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

sorted_gaia = gaia_results.copy()
sorted_gaia.sort('distance_ly_nominal', reverse=True)
top20 = sorted_gaia[:20]

###############################
# Part 2. Retrieve SIMBAD information for the top 20 objects.
###############################

# Set up SIMBAD
custom_simbad = Simbad()
custom_simbad.reset_votable_fields()
custom_simbad.add_votable_fields('MAIN_ID','ids','sp','otype')

def extract_luminosity_class(sp_type):
    if sp_type is None or sp_type.strip() == "":
        return "N/A"
    match = re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
    if match:
        return match.group(0)
    return "N/A"

def extract_common_name(ids_field, main_id):
    if ids_field is None:
        return "N/A"
    parts = [part.strip() for part in ids_field.split('|')]
    candidates = [x for x in parts if x != main_id]
    for candidate in candidates:
        if any(greek in candidate for greek in ['α','β','γ','δ','ε','ζ','η','θ','ι','κ','λ','μ','ν','ξ','ο','π','ρ','σ','τ','υ','φ','χ','ψ','ω']):
            return candidate
        if not re.match(r'^(HD|Gaia|TYC|2MASS|HIP)\s*\d+', candidate, re.IGNORECASE):
            return candidate
    return "N/A"

def get_simbad_info(ra, dec, radius=2*u.arcsec):
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    print(f"Querying SIMBAD at RA={ra:.5f}, Dec={dec:.5f}, radius={radius.to(u.arcsec)}")
    try:
        result = custom_simbad.query_region(coord, radius=radius)
    except Exception as e:
        print(f"SIMBAD query error at RA={ra:.5f}, Dec={dec:.5f}: {e}")
        result = None

    if result is None or len(result) == 0:
        if radius < 10*u.arcsec:
            return get_simbad_info(ra, dec, radius=radius + 3*u.arcsec)
        else:
            print(f"No SIMBAD match found for RA={ra:.5f}, Dec={dec:.5f} with radius up to {radius.to(u.arcsec)}")
            return {"main_id": "N/A", "common_name": "N/A", "sp_type": "N/A",
                    "lum_class": "N/A", "var_type": "N/A", "brightness_range": "N/A"}

    # Debug: print the columns available and the first row
    # print("SIMBAD result columns:", result.colnames)
    # print("SIMBAD first row:", result[0])
    
    main_id = result['MAIN_ID'][0] if 'MAIN_ID' in result.colnames else "N/A"
    if isinstance(main_id, bytes):
        main_id = main_id.decode('utf-8')
        
    ids_field = result['IDS'][0] if 'IDS' in result.colnames else "N/A"
    if isinstance(ids_field, bytes):
        ids_field = ids_field.decode('utf-8')
        
    sp_type = result['SP'][0] if 'SP' in result.colnames and result['SP'][0] is not None else "N/A"
    if isinstance(sp_type, bytes):
        sp_type = sp_type.decode('utf-8')
        
    lum_class = extract_luminosity_class(sp_type) if sp_type != "N/A" else "N/A"
    
    otype = result['OTYPE'][0] if 'OTYPE' in result.colnames and result['OTYPE'][0] is not None else "N/A"
    if isinstance(otype, bytes):
        otype = otype.decode('utf-8')
    var_type = otype if otype is not None and ("Var" in otype or "V*" in otype) else "N/A"
        
    brightness_range = "N/A"
    common_name = extract_common_name(ids_field, main_id)
    
    return {"main_id": main_id,
            "common_name": common_name,
            "sp_type": sp_type if sp_type is not None else "N/A",
            "lum_class": lum_class,
            "var_type": var_type,
            "brightness_range": brightness_range}

# (Optional) Test SIMBAD query with a known bright star (e.g., Sirius at RA=101.28716, Dec=-16.71612)
test_info = get_simbad_info(101.28716, -16.71612)
print("Test SIMBAD result for Sirius:", test_info)

sim_main_ids = []
sim_common_names = []
sim_sp_types = []
sim_lum_classes = []
sim_var_types = []
sim_brightness_ranges = []

print("\nQuerying SIMBAD for top-20 Gaia objects …")
for row in top20:
    ra = row['ra']
    dec = row['dec']
    sim_info = get_simbad_info(ra, dec)
    sim_main_ids.append(sim_info["main_id"])
    sim_common_names.append(sim_info["common_name"])
    sim_sp_types.append(sim_info["sp_type"])
    sim_lum_classes.append(sim_info["lum_class"])
    sim_var_types.append(sim_info["var_type"])
    sim_brightness_ranges.append(sim_info["brightness_range"])

top20['simbad_main_id'] = sim_main_ids
top20['common_name'] = sim_common_names
top20['spectral_type'] = sim_sp_types
top20['luminosity_class'] = sim_lum_classes
top20['variable_star_type'] = sim_var_types
top20['brightness_range'] = sim_brightness_ranges

###############################
# Part 3. Reorder and Format the Final Table
###############################

final_columns = ['common_name', 'simbad_main_id', 'source_id', 'phot_g_mean_mag',
                 'parallax', 'parallax_error',
                 'distance_ly_nominal', 'distance_ly_lower_bound',
                 'spectral_type', 'luminosity_class',
                 'variable_star_type', 'brightness_range']

final_table = top20[final_columns]

def format_row(row):
    return {
        'common_name': row['common_name'],
        'simbad_main_id': row['simbad_main_id'],
        'source_id': row['source_id'],
        'phot_g_mean_mag': f"{row['phot_g_mean_mag']:.3f}",
        'parallax': f"{row['parallax']:.5f}",
        'parallax_error': f"{row['parallax_error']:.5f}",
        'distance_ly_nominal': f"{row['distance_ly_nominal']:.5f}",
        'distance_ly_lower_bound': f"{row['distance_ly_lower_bound']:.5f}",
        'spectral_type': row['spectral_type'],
        'luminosity_class': row['luminosity_class'],
        'variable_star_type': row['variable_star_type'],
        'brightness_range': row['brightness_range']
    }

formatted_rows = [format_row(row) for row in final_table]
formatted_table = Table(rows=formatted_rows, names=final_columns)

print("\nFinal Table (top 20 Gaia objects with SIMBAD info):")
print(formatted_table)

formatted_table.write("gaia_top20_with_simbad.csv", format="csv", overwrite=True)
print("\nFinal table saved to 'gaia_top20_with_simbad.csv'.")
