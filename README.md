# URNG recurrences matrices

Linear recurrence proves to be a potent approach for uniform random number generation (URNG).  In 2005, David D. Thomas' publication titled "High Quality Uniform Random Number Generation Through LUT Optimized Linear Recurrences" illuminated an innovative path to optimize the utilization of lookup tables (LUTs) within such generators.

Based on Thomas' works, this project introduces a `search_matrix` function that efficiently locates an optimized linear recurrence matrix for any specified `k` (desired URN width) and `t` (number of LUT inputs).  Additionally, a table of such matrices `urng_rec_matrices.json` is provided, encompassing various frequently used combinations of `k` and `t` values.

In the JSON file, the first dimension is `k`, while the second dimension is `t`.  For each combination of `k` and `t` a matrix is stored in a sparse manner: It is stored as `k` lists, where each list is the indices of 1s in the `k`th row.
