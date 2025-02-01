#!/usr/bin/env python3
"""
Query Gaia DR3 for stars with Gaia G magnitude < 6, compute distances (in light years)
from the parallax and from (parallax + sigma), and print two sets of results:
1. Top 20 candidates based on nominal distance (naïve inversion of parallax)
2. Top 20 candidates based on a conservative lower–bound distance (using parallax + error)
"""

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.table import Table, vstack
import numpy as np

# Define the ADQL query.
# We select the source_id, phot_g_mean_mag, parallax and parallax_error.
# We restrict to sources with phot_g_mean_mag < 6 and parallax > 0.
query = """
SELECT source_id, phot_g_mean_mag, parallax, parallax_error
FROM gaiadr3.gaia_source
WHERE phot_g_mean_mag < 6
  AND parallax > 0
"""

print("Launching Gaia DR3 query …")
job = Gaia.launch_job(query)
results = job.get_results()
print("Query returned {} rows.".format(len(results)))

# Define conversion factor: 1 parsec = 3.26156 light years.
PC_TO_LY = 3.26156

def compute_distance_ly(parallax_mas):
    """Compute distance (in light years) using the naïve inversion: d = 1000 / parallax (mas)."""
    return (1000.0 / parallax_mas) * PC_TO_LY

def compute_lower_bound_distance_ly(parallax_mas, parallax_err_mas):
    """
    Compute a conservative lower bound on the distance (in light years) using the worst-case parallax:
    d_lb = 1000 / (parallax + error)
    """
    return (1000.0 / (parallax_mas + parallax_err_mas)) * PC_TO_LY

# Compute nominal distance and lower bound distance for each row.
nominal_distances = []
lower_bound_distances = []

for row in results:
    p = row['parallax']
    sigma = row['parallax_error']
    d_nom = compute_distance_ly(p)
    d_lb = compute_lower_bound_distance_ly(p, sigma)
    nominal_distances.append(d_nom)
    lower_bound_distances.append(d_lb)

# Add new columns to the results table.
results['distance_ly_nominal'] = nominal_distances
results['distance_ly_lower_bound'] = lower_bound_distances

# Create two sorted tables:
# 1. Sorted by nominal distance (largest first).
sorted_nominal = results.copy()
sorted_nominal.sort('distance_ly_nominal', reverse=True)
top20_nominal = sorted_nominal[:20]

# 2. Sorted by lower-bound distance (largest first).
sorted_lower_bound = results.copy()
sorted_lower_bound.sort('distance_ly_lower_bound', reverse=True)
top20_lower_bound = sorted_lower_bound[:20]

# Print the results in a formatted way.
def print_table(table, title):
    print("\n" + "="*len(title))
    print(title)
    print("="*len(title))
    # Print a header.
    header = ("source_id", "Gmag", "parallax (mas)", "σ (mas)", "d_nom (ly)", "d_lb (ly)")
    print("{:<25} {:>6} {:>15} {:>10} {:>12} {:>12}".format(*header))
    for row in table:
        print("{:<25} {:>6.2f} {:>15.3f} {:>10.3f} {:>12.0f} {:>12.0f}".format(
            str(row['source_id'])[:25],
            row['phot_g_mean_mag'],
            row['parallax'],
            row['parallax_error'],
            row['distance_ly_nominal'],
            row['distance_ly_lower_bound']
        ))

print_table(top20_nominal, "Top 20 Candidates by Nominal Distance")
print_table(top20_lower_bound, "Top 20 Candidates by Lower-Bound Distance (using parallax+σ)")

# Optionally, save the results to CSV files.
top20_nominal.write("gaia_top20_nominal.csv", format="csv", overwrite=True)
top20_lower_bound.write("gaia_top20_lower_bound.csv", format="csv", overwrite=True)
print("\nResults saved to 'gaia_top20_nominal.csv' and 'gaia_top20_lower_bound.csv'.")
