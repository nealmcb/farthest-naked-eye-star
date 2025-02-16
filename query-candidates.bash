# Run farthest_naked_eye_star.py on all the stars in candidates.csv
# Then join the results with candidates.csv to get the
# The latter step requires the Python csvkit:
#   https://csvkit.readthedocs.io/en/latest/scripts/csvlook.html.

# Pull the Gaia DR3 source_ids out from candidates.csv, and join with ','
ids="$(grep ,[0-9] candidates.csv | cut -d, -f2  | paste -sd ',' -))"

python farthest_naked_eye_star.py "source_id in ($ids)" \
 | tee candidates.log

mv gaia_top20_lowerbound.csv candidates-results.csv

csvjoin --left -c "Gaia DR3" candidates.csv candidates-results.csv > candidates-join.csv
