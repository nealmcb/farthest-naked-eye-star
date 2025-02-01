#!/usr/bin/env python3
"""
Query Gaia DR3 for stars with G < 6, compute distances (in light years)
from the parallax (and from parallax+error for a conservative lower bound),
then for the top 5 objects retrieve SIMBAD name and star–type info.
Final CSV headers are short:
  - common: Common name (if available)
  - main_id: SIMBAD main identifier
  - src_id: Gaia DR3 source_id
  - Gmag: Gaia phot_g_mean_mag
  - plx: parallax (mas)
  - plx_err: parallax_error (mas)
  - d_nom: nominal distance (ly) [integer]
  - d_lb: lower-bound distance (ly) [integer]
  - spec: spectral type
  - lum: luminosity class
  - var: variable star type
  - brange: brightness range

Distances are rounded to integer values.
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

gaia_query = """
SELECT source_id, ra, dec, phot_g_mean_mag, parallax, parallax_error
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

# Add distance columns (nominal and lower bound)
nominal_distances = []
lower_bound_distances = []
for row in gaia_results:
    p = row['parallax']
    sigma = row['parallax_error']
    nominal_distances.append(compute_distance_ly(p))
    lower_bound_distances.append(compute_lower_bound_distance_ly(p, sigma))

gaia_results['distance_ly_nominal'] = nominal_distances
gaia_results['distance_ly_lower_bound'] = lower_bound_distances

# Sort by nominal distance descending and take top 5 for debugging.
sorted_gaia = gaia_results.copy()
sorted_gaia.sort('distance_ly_nominal', reverse=True)
top5 = sorted_gaia[:5]

###############################
# Part 2. Retrieve SIMBAD Information for Top 5 Objects
###############################

custom_simbad = Simbad()
custom_simbad.reset_votable_fields()  # reset defaults
custom_simbad.add_votable_fields('MAIN_ID','ids','sp','otype')

def extract_luminosity_class(sp_type):
    if sp_type is None or sp_type.strip() == "":
        return "N/A"
    m = re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
    return m.group(0) if m else "N/A"

def extract_common_name(ids_field, main_id):
    if ids_field is None:
        return "N/A"
    parts = [part.strip() for part in ids_field.split('|')]
    # Remove main_id from list.
    candidates = [x for x in parts if x != main_id]
    for candidate in candidates:
        # If candidate contains Greek letters or isn't a catalog identifier.
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

    # Print the full first row of the result for debugging.
    print("SIMBAD result row:")
    for col in result.colnames:
        print(f"  {col}: {result[col][0]}")
    
    main_id = result['main_id'][0] if 'main_id' in result.colnames else "N/A"
    if isinstance(main_id, bytes):
        main_id = main_id.decode('utf-8')
    ids_field = result['ids'][0] if 'ids' in result.colnames else "N/A"
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
    brightness_range = "N/A"  # not provided by SIMBAD
    common_name = extract_common_name(ids_field, main_id)
    return {"main_id": main_id,
            "common_name": common_name,
            "sp_type": sp_type if sp_type is not None else "N/A",
            "lum_class": lum_class,
            "var_type": var_type,
            "brightness_range": brightness_range}

# Optional test: Query SIMBAD for Sirius (should return a match).
print("\nTesting SIMBAD query with Sirius (RA=101.28716, Dec=-16.71612):")
test_info = get_simbad_info(101.28716, -16.71612)
print("Test result for Sirius:", test_info)

# For each of the top 5 Gaia objects, query SIMBAD.
sim_main_ids = []
sim_common_names = []
sim_sp_types = []
sim_lum_classes = []
sim_var_types = []
sim_brightness_ranges = []

print("\nQuerying SIMBAD for top-5 Gaia objects …")
for row in top5:
    ra = row['ra']
    dec = row['dec']
    sim_info = get_simbad_info(ra, dec)
    sim_main_ids.append(sim_info["main_id"])
    sim_common_names.append(sim_info["common_name"])
    sim_sp_types.append(sim_info["sp_type"])
    sim_lum_classes.append(sim_info["lum_class"])
    sim_var_types.append(sim_info["var_type"])
    sim_brightness_ranges.append(sim_info["brightness_range"])

# For debugging, print out the SIMBAD info lists.
print("\nSIMBAD main IDs for top 5:", sim_main_ids)
print("SIMBAD common names for top 5:", sim_common_names)
print("SIMBAD spectral types for top 5:", sim_sp_types)

###############################
# Part 3. Merge Results into a Final Table with Shortened Headers
###############################

# We'll construct a new list of dictionaries (one per object) with these headers:
# common, main_id, src_id, Gmag, plx, plx_err, d_nom, d_lb, spec, lum, var, brange
final_data = []
for i, row in enumerate(top5):
    # Round distances to integers.
    d_nom_int = int(round(row['distance_ly_nominal']))
    d_lb_int = int(round(row['distance_ly_lower_bound']))
    final_data.append({
        'common': sim_common_names[i],
        'main_id': sim_main_ids[i],
        'src_id': row['source_id'],
        'Gmag': f"{row['phot_g_mean_mag']:.3f}",
        'plx': f"{row['parallax']:.5f}",
        'plx_err': f"{row['parallax_error']:.5f}",
        'd_nom': f"{d_nom_int}",
        'd_lb': f"{d_lb_int}",
        'spec': sim_sp_types[i],
        'lum': sim_lum_classes[i],
        'var': sim_var_types[i],
        'brange': sim_brightness_ranges[i]
    })

final_table = Table(rows=final_data, names=['common','main_id','src_id','Gmag','plx','plx_err','d_nom','d_lb','spec','lum','var','brange'])

print("\nFinal Table (top 5 Gaia objects with SIMBAD info):")
print(final_table)

# Save the final table to a CSV file.
final_table.write("gaia_top5_with_simbad.csv", format="csv", overwrite=True)
print("\nFinal table saved to 'gaia_top5_with_simbad.csv'.")
