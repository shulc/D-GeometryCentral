#ifndef GC_C_H
#define GC_C_H
#ifdef __cplusplus
extern "C" {
#endif

// Heat method geodesic distance
// verts: flat xyz positions, nv vertices
// tris: flat triangle indices (3 per tri), nt triangles
// sources: source vertex indices, ns sources
// out_dist: output distances, size nv
// Returns 0 on success, non-zero on error
int gc_geodesic_distance(const double* verts, int nv, const unsigned int* tris, int nt,
                         const int* sources, int ns, double* out_dist);

// Disk parameterization via BFF (Boundary First Flattening — conformal, free boundary)
// Requires a mesh with at least one boundary loop.
// verts: flat xyz positions, nv vertices
// tris: flat triangle indices, nt triangles
// out_uv: output UV coordinates, size 2*nv (u0,v0, u1,v1, ...)
// Returns 0 on success, 1 on internal error, 3 if mesh has no boundary
int gc_parameterize_disk(const double* verts, int nv, const unsigned int* tris, int nt,
                         double* out_uv);

// Cotan-Laplacian smoothing (boundary pinned)
// verts: flat xyz positions, nv vertices
// tris: flat triangle indices, nt triangles
// iters: number of smoothing iterations
// lambda: step size
// out_verts: output positions, size 3*nv
// Returns 0 on success, non-zero on error
int gc_cotan_smooth(const double* verts, int nv, const unsigned int* tris, int nt,
                    int iters, double lambda, double* out_verts);

#ifdef __cplusplus
}
#endif
#endif
