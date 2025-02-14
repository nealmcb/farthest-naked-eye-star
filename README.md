# Furthest naked-eye stars via Gaia

Python software to search Gaia for "most distant naked-eye star" candidates

Code written mostly by ChatGPT o3-mini, evolved via dozens of prompts by Neal McBurnett, 2025-02-01

# Usage
```
python furthest_naked_eye_star.py  | tee output-log.txt
csvlook -I gaia_top20_lowerbound.csv
```

For line execution and coverage counts:

```
python -m trace --count --coverdir=trace_results furthest_naked_eye_star.py
more trace_results/furthest_visible_star.cover
```

# TODO
* Can SIMBAD query be sped up or run in parallel via Gaia IDs for these rather bright stars?
  Individual queries take perhaps 10 seconds
* Confirm that get_simbad_info() is getting the right stars
* Can extract_common_name() be refined?
* Clean up code by moving functions to the top
* Make it usable as a library by moving code to main()
* Validate spectral type filtering: re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
  But note that so far they all passed

* Provide answers for a given location for a given time and weather conditions
* Revisit when Gaia DR4 comes out, with significantly improved distances
