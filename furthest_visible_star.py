#!/usr/bin/env python3
"""
This script queries Gaia DR3 for stars with G < 6 and positive parallax,
computes two distance estimates (nominal and conservative lower-bound) in light years,
and then produces two sets of results:
  - top 20 by nominal (expected) distance
  - top 20 by lower-bound distance
For each candidate the script queries SIMBAD to obtain:
  - SIMBAD main identifier (main_id)
  - A common name (ignoring Gaia names)
  - Visual magnitude (Vmag, from flux(V))
  - Spectral type (sp_type) and a parsed luminosity class (lum)
  - Variable star type (var)
  - (Brightness range is set to "N/A")
Finally, the output table (with shortened headers) is constructed with the following columns:
  1. main_id (SIMBAD main identifier)
  2. common (common name)
  3. Gaia DR3 (Gaia source id)
  4. Vmag (SIMBAD visual magnitude)
  5. Gmag (Gaia phot_g_mean_mag)
  6. ly_nom (nominal distance in ly, rounded to integer)
  7. ly_lb (lower-bound distance in ly, rounded to integer)
  8. spec (spectral type)
  9. lum (luminosity class)
 10. var (variable star type)
 11. brange (brightness range, "N/A")
 12. Con (constellation abbreviation, computed from RA,Dec)
 13. RA (Gaia right ascension)
 14. Dec (Gaia declination)
Two CSV files are saved:
  - gaia_top20_nominal.csv
  - gaia_top20_lowerbound.csv
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

# Query Gaia DR3 (with a TOP clause to get up to 10,000 rows)
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

# Conversion factor: 1 pc = 3.26156 ly.
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

# Create two sorted copies:
sorted_nom = gaia_results.copy()
sorted_nom.sort('distance_ly_nominal', reverse=True)
top20_nom = sorted_nom[:20]

sorted_lb = gaia_results.copy()
sorted_lb.sort('distance_ly_lower_bound', reverse=True)
top20_lb = sorted_lb[:20]

###############################
# Part 2. SIMBAD Query Function and Helpers
###############################

# Set up SIMBAD: request MAIN_ID, ids, sp_type, otype, and flux(V) for visual mag.
custom_simbad = Simbad()
custom_simbad.reset_votable_fields()
custom_simbad.add_votable_fields('MAIN_ID','ids','sp_type','otype','flux(V)')

def extract_luminosity_class(sp_type):
    if sp_type is None or sp_type.strip() == "":
        return "N/A"
    m = re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
    return m.group(0) if m else "N/A"

def extract_common_name(ids_field, main_id):
    if ids_field is None:
        return "N/A"
    parts = [part.strip() for part in ids_field.split('|')]
    # Exclude the main_id and any candidate containing "Gaia"
    candidates = [x for x in parts if x != main_id and "gaia" not in x.lower()]
    for candidate in candidates:
        if any(greek in candidate for greek in ['α','β','γ','δ','ε','ζ','η','θ','ι','κ','λ','μ','ν','ξ','ο','π','ρ','σ','τ','υ','φ','χ','ψ','ω']):
            return candidate
        if not re.match(r'^(HD|TYC|2MASS|HIP)\s*\d+', candidate, re.IGNORECASE):
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

    # Normalize column names to lowercase using rename_column.
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
                    "lum_class": "N/A", "var_type": "N/A", "brightness_range": "N/A", "vmag": "N/A"}

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
    
    # Get visual magnitude: check for column "v" instead of "flux(v)"
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
# Part 3. Enrich Each Top 20 Set with SIMBAD Data
###############################

def enrich_with_simbad(top_table):
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

print("\nEnriching top20 nominal candidates with SIMBAD info …")
top20_nom = enrich_with_simbad(top20_nom)
print("\nEnriching top20 lower-bound candidates with SIMBAD info …")
top20_lb = enrich_with_simbad(top20_lb)

###############################
# Part 4. Build Final Tables with New Column Order and Labels
###############################

# We'll construct the final table with these columns:
# 1. main_id, 2. common, 3. Gaia DR3, 4. Vmag, 5. Gmag, 6. ly_nom, 7. ly_lb,
# 8. spec, 9. lum, 10. var, 11. brange, 12. Con, 13. RA, 14. Dec
def build_final_table(top_table):
    final_data = []
    for row in top_table:
        d_nom_int = int(round(row['distance_ly_nominal']))
        d_lb_int = int(round(row['distance_ly_lower_bound']))
        # Compute constellation abbreviation from RA, Dec.
        coord = SkyCoord(ra=row['ra']*u.deg, dec=row['dec']*u.deg, frame='icrs')
        constellation = coord.get_constellation(short_name=True)
        final_data.append({
            'main_id': row['simbad_main_id'],         # SIMBAD main identifier
            'common': row['common_name'],             # common name
            'Gaia DR3': row['source_id'],             # Gaia source id
            'Vmag': f"{float(row['vmag']):.3f}" if row['vmag']!="N/A" else "N/A",
            'Gmag': f"{row['phot_g_mean_mag']:.3f}",
            'ly_nom': f"{d_nom_int}",
            'ly_lb': f"{d_lb_int}",
            'spec': row['sp_type'],
            'lum': row['lum_class'],
            'var': row['var_type'],
            'brange': row['brightness_range'],
            'Con': constellation,
            'RA': f"{row['ra']:.5f}",
            'Dec': f"{row['dec']:.5f}"
        })
    # Define final headers (order as desired)
    return Table(rows=final_data, names=['main_id','common','Gaia DR3','Vmag','Gmag','ly_nom','ly_lb','spec','lum','var','brange','Con','RA','Dec'])

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
