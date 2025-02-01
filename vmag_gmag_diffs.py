import pandas as pd

# Change the filename to your CSV file as needed.
filename = "gaia_top20_nominal.csv"

# Read the CSV file.
df = pd.read_csv(filename)

# Remove rows where Vmag is "N/A" or Gmag is "N/A" and convert columns to floats.
df = df[(df["Vmag"] != "N/A") & (df["Gmag"] != "N/A")].copy()
df["Vmag"] = pd.to_numeric(df["Vmag"], errors='coerce')
df["Gmag"] = pd.to_numeric(df["Gmag"], errors='coerce')

# Compute the difference (Visual magnitude minus Gaia G magnitude).
df["diff"] = df["Vmag"] - df["Gmag"]

# Print descriptive statistics of the difference.
print("Descriptive statistics for Vmag - Gmag:")
print(df["diff"].describe())
