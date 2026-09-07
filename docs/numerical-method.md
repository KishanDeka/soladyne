# Translation map

`Model::build_profiles()` corresponds to Fortran Parts I, II, and II-A.
`Model::load_initial_state()` corresponds to Part III. `build_coefficients()`
is Part IV. `advance_one()` retains the statement-group order of Part V-A;
`apply_boundaries()` is Part V-B. The active output blocks are implemented by
`step_with_output()` and `write_final_outputs()`.

One deterministic clarification was required: the Fortran evaluates several
midpoint helper expressions at `sin(theta)=0`, where uninitialized boundary
array values can participate in `0/0`. Those endpoint helper entries are set to
zero in C++; the ADI stencil only consumes their analytically finite products.
No interior formula, solve, or active physical term was removed.

The C++ grids deliberately use 1-based indexing internally. This makes every
Fortran index—including `2*i-1`, `2*j-3`, and `n+1`—directly auditable against
the source and avoids error-prone index algebra during translation. Snapshot
arrays returned to Python are conventional contiguous `(257,257)` arrays in
the same radial-major, theta-minor order as the text files.
