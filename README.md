# D-GeometryCentral

D bindings for a thin C shim around **geometry-central** — a C++ library for
computational geometry on surface meshes by Nicholas Sharp et al.

Three operations are exposed:

| Entry point | What it does |
|---|---|
| `gc_geodesic_distance` | Heat-method geodesic distance from a set of source vertices |
| `gc_parameterize_disk` | Conformal disk parameterization (UV flattening, no cuts required) |
| `gc_cotan_smooth` | Explicit cotan-Laplacian smoothing (boundary vertices pinned) |

## Layout

```
include/gc_c.h                 — C API (extern "C")
csrc/gc_c.cpp                  — implementation over geometry-central
source/geometrycentral/c.d     — D extern(C) declarations, 1:1 with gc_c.h
CMakeLists.txt                 — builds libgc_c.a + vendored geometry-central
extern/geometry-central        — git submodule, pinned to commit 019669d (v1.1.0 + 3 commits)
examples/geodesic.d            — smoke test: geodesics + parameterize + smooth
```

## First build

```sh
git submodule update --init --recursive
dub build
```

The `preBuildCommands-posix` hook runs CMake on the submodule + shim before
dub compiles the D side. geometry-central fetches Eigen and Spectra via
CMake FetchContent on first configure (~1–2 min on first run, cached after).
The full first cmake build takes ~5–15 min depending on machine.

## Smoke test

```sh
dub run --config=geodesic
```

Expected output (exact numbers vary slightly by platform):

```
Icosahedron: 12 vertices, 20 triangles

[Geodesic distance from vertex 0]  rc=0
  dist[0] = 0.000000  (expected 0.0)
  dist[1] = 1.…
  max dist = 3.…  (sphere π ≈ 3.141593)
  OK

[Geodesic distance from vertices 0+3]  rc=0
  dist[0] = 0.000000  dist[3] = 0.000000
  OK

[Parameterize disk — flat quad grid 3x3 = 4 quads = 8 tris]
  grid: 9 verts, 8 tris
  rc = 0
  …
  OK

[Cotan smooth — 3x3 grid, 5 iters, lambda=0.5]
  rc = 0
  centre z before=1.0000  after=…  (smoothed toward 0)
  corner[0] z before=0.0000  after=0.0000  (boundary pinned)
  OK

All tests PASSED
```

## Consuming from another dub project

```json
"dependencies": {
    "d-geometrycentral": { "path": "../D-GeometryCentral" }
}
```

Then `import geometrycentral.c;`.

## API

```d
import geometrycentral.c;

// Geodesic distance
auto dist = new double[](nv);
int[1] src = [0];
gc_geodesic_distance(verts.ptr, nv, tris.ptr, nt, src.ptr, 1, dist.ptr);

// Disk parameterization
auto uv = new double[](2 * nv);
gc_parameterize_disk(verts.ptr, nv, tris.ptr, nt, uv.ptr);

// Cotan smoothing
auto smoothed = new double[](3 * nv);
gc_cotan_smooth(verts.ptr, nv, tris.ptr, nt, 10, 0.5, smoothed.ptr);
```

## License

The shim + bindings are MIT-licensed. geometry-central upstream is also MIT.
Eigen and Spectra (fetched at build time) are MPL-2.0.
