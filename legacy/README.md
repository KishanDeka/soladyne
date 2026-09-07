# Legacy sources

`source/surya_dynamo_v3.f` is the authoritative source for the exact C++ port.
The earlier v2 and omega-only variants plus their plotting scripts are retained
for provenance, but they do not define the v3 implementation.

The C++ translation keeps active v3 statements and file contracts intact while
moving common-block/global state into `solar_dynamo::Model`. Commented diagnostics
and buoyancy experiments remain inactive, as in the supplied v3 source.
