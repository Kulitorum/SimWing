# AMReX provenance

SimWing's optional production-grid evaluation uses AMReX 26.06 as an external
CMake dependency. The default FetchContent path is pinned to upstream commit
`bf15fdce52093c715a1dd9a9da64ba66f0233be1`; users may instead provide a
compatible installed package through `SIMWING_AMREX_ROOT`.

Upstream: <https://github.com/AMReX-Codes/amrex>

No AMReX source is vendored here. `LICENSE` and `NOTICE` reproduce the
upstream attribution shipped with the pinned dependency for source and binary
redistribution.
