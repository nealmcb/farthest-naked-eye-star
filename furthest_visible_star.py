#!/usr/bin/env python3
"""
Query Gaia DR3 for stars with G < 6, compute distances (in light years)
from the parallax (and from parallax+error for a conservative lower bound),
then for the top 20 objects (by distance) retrieve SIMBAD name and star–type info.
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

# ADQL query: retrieve Gaia DR3 sources with G < 6 and positive parallax.
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

# Conversion: 1 parsec = 3.26156 light years.
PC_TO_LY = 3.26156

def compute_distance_ly(parallax_mas):
    """Naïve distance (in light years) from parallax (in mas)."""
    return (1000.0 / parallax_mas) * PC_TO_LY

def compute_lower_bound_distance_ly(parallax_mas, parallax_err_mas):
    """Conservative lower-bound distance (in ly) using (parallax + error)."""
    return (1000.0 / (parallax_mas + parallax_err_mas)) * PC_TO_LY

# Compute distances for each star and add as new columns.
nominal_distances = []
lower_bound_distances = []
for row in gaia_results:
    p = row['parallax']
    sigma = row['parallax_error']
    d_nom = compute_distance_ly(p)
    d_lb = compute_lower_bound_distance_ly(p, sigma)
    nominal_distances.append(d_nom)
    lower_bound_distances.append(d_lb)

gaia_results['distance_ly_nominal'] = nominal_distances
gaia_results['distance_ly_lower_bound'] = lower_bound_distances

# For our purposes, we choose the top 20 objects sorted by nominal distance (largest first).
sorted_gaia = gaia_results.copy()
sorted_gaia.sort('distance_ly_nominal', reverse=True)
top20 = sorted_gaia[:20]

###############################
# Part 2. Retrieve SIMBAD information for the top 20 objects.
###############################

# Set up a custom SIMBAD query that returns:
# - MAIN_ID (SIMBAD main identifier)
# - IDS (all alternative identifiers)
# - SP_TYPE (spectral type)
# - OTYPE (object type, which we will use to help flag variable stars)
custom_simbad = Simbad()
custom_simbad.remove_votable_fields('coordinates')  # already have RA,Dec
custom_simbad.add_votable_fields('ids', 'sp', 'otype')

def extract_luminosity_class(sp_type):
    """
    Try to extract the luminosity class from a spectral type string.
    This function looks for common luminosity classes (I, II, III, IV, V) possibly with additional letters.
    If none is found, return "N/A".
    """
    if sp_type is None or sp_type.strip() == "":
        return "N/A"
    # Example: "M2Iab" -> look for I, II, III, IV, or V patterns
    match = re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
    if match:
        return match.group(0)
    return "N/A"

def extract_common_name(ids_field, main_id):
    """
    Given the SIMBAD IDS field (a string of alternative names, separated by |)
    and the main_id, try to pick one common name that isn't just a catalog number.
    This is heuristic; here we choose the first identifier that doesn't start with "Gaia" or "HD" or similar.
    If no such candidate is found, return "N/A".
    """
    if ids_field is None:
        return "N/A"
    # Split the IDS string (identifiers are separated by '|' in SIMBAD)
    parts = [part.strip() for part in ids_field.split('|')]
    # Remove the main_id from the list
    candidates = [x for x in parts if x != main_id]
    # Look for a candidate that looks like a common name:
    for candidate in candidates:
        # Heuristic: if it contains a Greek letter or is not just a catalog number.
        # For example, "Alpha Centauri" or "Sirius" are good.
        if any(greek in candidate for greek in ['α','β','γ','δ','ε','ζ','η','θ','ι','κ','λ','μ','ν','ξ','ο','π','ρ','σ','τ','υ','φ','χ','ψ','ω']):
            return candidate
        # Also, if the candidate contains a common name-like string (e.g., "Capella")
        if not re.match(r'^(HD|Gaia|TYC|2MASS|HIP)\s*\d+', candidate, re.IGNORECASE):
            return candidate
    return "N/A"

def get_simbad_info(ra, dec):
    """
    Given coordinates (in degrees), do a cone search in SIMBAD (radius 2 arcsec)
    and return a dict with the following keys:
       - main_id
       - common_name
       - sp_type
       - lum_class
       - var_type
       - brightness_range
    If no result is found, return "N/A" for all.
    """
    coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    try:
        result = custom_simbad.query_region(coord, radius=2*u.arcsec)
    except Exception as e:
        result = None

    if result is None or len(result) == 0:
        return {"main_id": "N/A", "common_name": "N/A", "sp_type": "N/A",
                "lum_class": "N/A", "var_type": "N/A", "brightness_range": "N/A"}

    # Use the first match.
    main_id = result['MAIN_ID'][0].decode('utf-8') if isinstance(result['MAIN_ID'][0], bytes) else result['MAIN_ID'][0]
    ids_field = result['IDS'][0].decode('utf-8') if isinstance(result['IDS'][0], bytes) else result['IDS'][0]
    sp_type = result['SP_TYPE'][0].decode('utf-8') if (result['SP_TYPE'][0] is not None and isinstance(result['SP_TYPE'][0], bytes)) else (result['SP_TYPE'][0] if result['SP_TYPE'][0] is not None else "N/A")
    lum_class = extract_luminosity_class(sp_type) if sp_type != "N/A" else "N/A"
    # For variable star type, we check the object type field (OTYPE).
    otype = result['OTYPE'][0].decode('utf-8') if isinstance(result['OTYPE'][0], bytes) else result['OTYPE'][0]
    if otype is not None and ("Var" in otype or "V*" in otype):
        var_type = otype
    else:
        var_type = "N/A"
    # Brightness range is not typically available in SIMBAD; we fill "N/A"
    brightness_range = "N/A"
    common_name = extract_common_name(ids_field, main_id)
    return {"main_id": main_id,
            "common_name": common_name,
            "sp_type": sp_type if sp_type is not None else "N/A",
            "lum_class": lum_class,
            "var_type": var_type,
            "brightness_range": brightness_range}

# For the top20 Gaia objects, query SIMBAD and build lists for the new columns.
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

# Add these SIMBAD fields as new columns to the top20 table.
top20['simbad_main_id'] = sim_main_ids
top20['common_name'] = sim_common_names
top20['spectral_type'] = sim_sp_types
top20['luminosity_class'] = sim_lum_classes
top20['variable_star_type'] = sim_var_types
top20['brightness_range'] = sim_brightness_ranges

###############################
# Part 3. Reorder and Output the Final Table
###############################

# Desired column order:
#   1. common_name
#   2. simbad_main_id
#   3. Gaia source_id
#   4. phot_g_mean_mag
#   5. parallax
#   6. parallax_error
#   7. distance_ly_nominal
#   8. distance_ly_lower_bound
#   9. spectral_type
#  10. luminosity_class
#  11. variable_star_type
#  12. brightness_range

final_columns = ['common_name', 'simbad_main_id', 'source_id', 'phot_g_mean_mag',
                 'parallax', 'parallax_error',
                 'distance_ly_nominal', 'distance_ly_lower_bound',
                 'spectral_type', 'luminosity_class',
                 'variable_star_type', 'brightness_range']

final_table = top20[final_columns]

# Print the final table with a formatted header.
print("\nFinal Table (top 20 Gaia objects with SIMBAD info):")
print(final_table)

# Optionally, save the final table to a CSV file.
final_table.write("gaia_top20_with_simbad.csv", format="csv", overwrite=True)
print("\nFinal table saved to 'gaia_top20_with_simbad.csv'.")
