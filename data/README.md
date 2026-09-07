# Legacy runtime data

The exact solver reads one file: `init.dat`, containing 66,049 rows and four
Fortran-formatted columns `(theta, radius, r*sin(theta)*u, ub)`. Put it here or
pass its path to `Model.initialize()` / `solar-dynamo --init`.

The supplied files were intentionally not duplicated in the source archive because
`psi(1).dat` is 64 MB and contains six identical blocks. Recorded SHA-256 checksums:

| supplied file | SHA-256 |
|---|---|
| `init(1).dat` | `26e21ab08cf13b24a028979004be575ffaa41b102d41a780f093601c4847c02a` |
| `psi(1).dat` | `1404a76f85db55cad7f08c8880923a62f433bba1f9b2bb693e696cb893736725` |
| `diffrot(1).dat` | `403e8f4859ba23ad0c4c4d019babf222b792de39c09c643295cb9ff160ff1f78` |
| `omg_lat(1).dat` | `a750dcdcd26ce50edaf993e02bd57d864f9656dfd3a1b4520a5cc693f8e3919e` |

`psi.dat` is not a runtime input: v3 calculates the stream function internally.
The supplied `omg_lat` contains three appended run segments (442, 95, and 474
snapshots); the last two use 0.008 output spacing and therefore do not match the
unmodified v3 `dt=0.00005` configuration.
