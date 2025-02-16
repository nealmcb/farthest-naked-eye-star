# Farthest naked-eye stars via Gaia

Python software to search Gaia for "most distant naked-eye star" candidates

Code written with significant help from ChatGPT o3-mini, evolved and polished
after dozens of prompts by Neal McBurnett, 2025-02

# Usage
```
python farthest_naked_eye_star.py  | tee output-log.txt
csvlook -I gaia_top20_lowerbound.csv
```

csvlook is part of [csvkit](https://csvkit.readthedocs.io/en/latest/index.html)

The code also allows for custom queries, like this:

`python farthest_naked_eye_star.py "source_id = 2005992002061917312"`

and as demonstrated in the `query-candidates.bash` script.

# Notes
Thanks to the contributors to the discussions at 
[What is the farthest-away star visible to the naked eye? Physics Stackexchange](https://physics.stackexchange.com/questions/45759/what-is-the-farthest-away-star-visible-to-the-naked-eye)
and Cloudy Nights,
[What is the farthest star visible naked\-eye? \- Deep Sky Observing \- Cloudy Nights](https://www.cloudynights.com/topic/623558-what-is-the-farthest-star-visible-naked-eye/) (free login needed),
for valuable insights and suggesting a variety of good candidate stars.

# Debugging
For line execution and coverage counts:

```
python -m trace --count --coverdir=trace_results farthest_naked_eye_star.py
more trace_results/farthest_visible_star.cover
```

# TODO
* Can SIMBAD query be sped up or run in parallel via Gaia IDs for these rather bright stars?
  Individual queries take perhaps 10 seconds
* Confirm that get_simbad_info() is getting the right stars
* Can extract_common_name() be refined?
* Clean up code by moving functions to the top
* Make it usable as a library by moving code to main()
* For custom queries, either sort while including unsortable rows also, or leave unsorted but only process one table
* Validate spectral type filtering: re.search(r'(I{1,3}[ab]?)|(IV)|(V)', sp_type)
  But note that so far they all passed

* Provide answers for a given location for a given time and weather conditions
* Revisit when Gaia DR4 comes out, with significantly improved distances
