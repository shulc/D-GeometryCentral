/// D bindings for the D-GeometryCentral C shim.
///
/// Mirrors include/gc_c.h exactly — one extern(C) declaration per
/// entry point. Higher-level D wrappers (RAII, slice-based API, GC-safe
/// wrappers) belong in sibling modules.
module geometrycentral.c;

extern(C) @nogc nothrow:

/// Compute heat-method geodesic distance from a set of source vertices.
///
/// Params:
///   verts    = flat xyz positions [x0,y0,z0, x1,...], length = 3*nv
///   nv       = number of vertices
///   tris     = flat triangle indices [i0,i1,i2, ...], length = 3*nt
///   nt       = number of triangles
///   sources  = source vertex indices, length = ns
///   ns       = number of source vertices
///   out_dist = output distances per vertex, caller-allocated length nv
/// Returns 0 on success, 1 on exception, 2 on out-of-range source index.
int gc_geodesic_distance(const(double)* verts, int nv,
                          const(uint)*   tris,  int nt,
                          const(int)*    sources, int ns,
                          double* out_dist);

/// Flatten a mesh with boundary to 2D UV coordinates via BFF
/// (Boundary First Flattening — conformal, free boundary).
/// Requires at least one boundary loop.
///
/// Params:
///   verts  = flat xyz positions, length = 3*nv
///   nv     = number of vertices
///   tris   = flat triangle indices, length = 3*nt
///   nt     = number of triangles
///   out_uv = output UV coords [u0,v0, u1,v1, ...], caller-allocated length 2*nv
/// Returns 0 on success, 1 on internal error, 3 if mesh has no boundary.
int gc_parameterize_disk(const(double)* verts, int nv,
                          const(uint)*   tris,  int nt,
                          double* out_uv);

/// Apply explicit cotan-Laplacian smoothing (boundary vertices pinned).
///
/// Params:
///   verts     = flat xyz positions, length = 3*nv
///   nv        = number of vertices
///   tris      = flat triangle indices, length = 3*nt
///   nt        = number of triangles
///   iters     = number of smoothing iterations
///   lambda    = step size (0 < lambda <= 1 for stability)
///   out_verts = output positions, caller-allocated length 3*nv
/// Returns 0 on success, non-zero on error.
int gc_cotan_smooth(const(double)* verts, int nv,
                     const(uint)*   tris,  int nt,
                     int iters, double lambda,
                     double* out_verts);
