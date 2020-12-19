# Nussinov Algorithm: Project for UIUC CS466 Introduction to Bioinformatics

For a RNA sequence, the Nussinov algorithm can find the pseudoknot-free secondary structure with maximum number of complementary base  pairings in O(n^3) time, where n is the length of the input RNA sequence.

This repo provides the Python implementation to:
1. the Nussinov algorithm to fill the DP table to find the maximum number of pseudoknot-free complementary base pairs.
2. reconstruct the secondary structure by tracing back from the constructed table.
3. visualize the backtrace according to the DP table and reconstruction.

## Dependencies
* matplotlib=3.3.3
* numpy=1.19.2

## Running Visualization on Single Sequence
First, change the configurations in `vis_config.py`, where
* `pairing`: a list of tuples specifying complementary base pairs allowed
* `hairpin`: int, minimum length of hairpin loops
* `max_num_sol`: int, maximum number of solutions returned by the reconstruction algorithm (or set to `None` for all solutions)
* `output_dir`: string, where the visualization output of the reconstruction is saved

Then, run
```
python visualize_main.py -- input <INPUT_STRING>
```

The command generates one plot of the backtrace for each reconstruction of the input string. We provide the output for the sequence `GGUCCAC` with allowed pairings of A-U, C-G, and G-U, and the minimum length of hairpin loop is set to 1 according to the example in [the Freiburg RNA Tool](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Nussinov).

The algorithm finds 5 solutions with a score of 2:
![traceback_((.))...png](./output/GGUCCAC_A-U_C-G_G-U_hairpin_1/traceback_((.))...png)

## Running Test Suites
TODO

## Authors
- [Zhepei Wang](https://github.com/zhepeiw)
- [Vivian Hu](https://github.com/vivianh2)
