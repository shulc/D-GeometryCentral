#include "gc_c.h"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/boundary_first_flattening.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <memory>
#include <vector>
#include <stdexcept>
#include <cstdio>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ---------------------------------------------------------------------------
// Helper: build ManifoldSurfaceMesh + VertexPositionGeometry from flat arrays
// ---------------------------------------------------------------------------
static std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
                 std::unique_ptr<VertexPositionGeometry>>
buildMesh(const double* verts, int nv, const unsigned int* tris, int nt)
{
    std::vector<std::vector<size_t>> faces(nt);
    for (int i = 0; i < nt; i++)
        faces[i] = {(size_t)tris[3*i], (size_t)tris[3*i+1], (size_t)tris[3*i+2]};

    auto mesh = std::make_unique<ManifoldSurfaceMesh>(faces);

    VertexData<Vector3> positions(*mesh);
    for (int i = 0; i < nv; i++)
        positions[mesh->vertex(i)] = Vector3{verts[3*i], verts[3*i+1], verts[3*i+2]};

    auto geom = std::make_unique<VertexPositionGeometry>(*mesh, positions);
    return {std::move(mesh), std::move(geom)};
}

// ---------------------------------------------------------------------------
// Heat-method geodesic distance from a set of source vertices
// ---------------------------------------------------------------------------
extern "C" int gc_geodesic_distance(const double* verts, int nv, const unsigned int* tris, int nt,
                                     const int* sources, int ns, double* out_dist)
{
    try {
        auto [mesh, geom] = buildMesh(verts, nv, tris, nt);

        HeatMethodDistanceSolver solver(*geom);

        std::vector<Vertex> srcVerts;
        srcVerts.reserve(ns);
        for (int i = 0; i < ns; i++) {
            if (sources[i] < 0 || sources[i] >= nv)
                return 2; // out-of-range source
            srcVerts.push_back(mesh->vertex(sources[i]));
        }

        VertexData<double> dist = solver.computeDistance(srcVerts);
        for (int i = 0; i < nv; i++)
            out_dist[i] = dist[mesh->vertex(i)];
        return 0;
    } catch (...) { return 1; }
}

// ---------------------------------------------------------------------------
// Disk parameterization — BFF (Boundary First Flattening)
// Works on any mesh with at least one boundary loop. The boundary is
// flattened to a free-boundary conformal map (not constrained to a circle).
// ---------------------------------------------------------------------------
extern "C" int gc_parameterize_disk(const double* verts, int nv, const unsigned int* tris, int nt,
                                     double* out_uv)
{
    try {
        auto [mesh, geom] = buildMesh(verts, nv, tris, nt);

        if (!mesh->hasBoundary())
            return 3; // BFF requires a mesh with boundary

        VertexData<Vector2> uvs = parameterizeBFF(*mesh, *geom);
        for (int i = 0; i < nv; i++) {
            Vector2 uv = uvs[mesh->vertex(i)];
            out_uv[2*i]   = uv.x;
            out_uv[2*i+1] = uv.y;
        }
        return 0;
    } catch (const std::exception& e) {
        fprintf(stderr, "[gc_parameterize_disk] exception: %s\n", e.what());
        return 1;
    } catch (...) { return 1; }
}

// ---------------------------------------------------------------------------
// Explicit cotan-Laplacian smoothing (boundary pinned via skip-boundary logic)
//
// Each iteration: new_pos[v] = pos[v] + lambda * (weighted_avg_neighbors - pos[v])
// where weights are cotan weights, boundary vertices kept fixed.
// ---------------------------------------------------------------------------
extern "C" int gc_cotan_smooth(const double* verts, int nv, const unsigned int* tris, int nt,
                                int iters, double lambda, double* out_verts)
{
    try {
        auto [mesh, geom] = buildMesh(verts, nv, tris, nt);

        // Require halfedge cotan weights
        geom->requireHalfedgeCotanWeights();
        const HalfedgeData<double>& cotanW = geom->halfedgeCotanWeights;

        // Work in a mutable positions array
        std::vector<Vector3> pos(nv);
        for (int i = 0; i < nv; i++)
            pos[i] = Vector3{verts[3*i], verts[3*i+1], verts[3*i+2]};

        for (int iter = 0; iter < iters; iter++) {
            std::vector<Vector3> newPos = pos;
            for (Vertex v : mesh->vertices()) {
                if (v.isBoundary()) continue; // pin boundary

                double wSum = 0.0;
                Vector3 laplacian{0, 0, 0};
                for (Halfedge he : v.outgoingHalfedges()) {
                    // cotan weight for the twin halfedge (twin's cotan is
                    // the weight associated with the edge from v's perspective)
                    double w = cotanW[he] + cotanW[he.twin()];
                    if (w < 0) w = 0; // clamp negative weights (non-Delaunay)
                    Vertex nb = he.tipVertex();
                    laplacian += w * (pos[nb.getIndex()] - pos[v.getIndex()]);
                    wSum += w;
                }
                if (wSum > 1e-12)
                    newPos[v.getIndex()] = pos[v.getIndex()] + lambda * (laplacian / wSum);
            }
            pos = newPos;
        }

        for (int i = 0; i < nv; i++) {
            out_verts[3*i]   = pos[i].x;
            out_verts[3*i+1] = pos[i].y;
            out_verts[3*i+2] = pos[i].z;
        }
        return 0;
    } catch (...) { return 1; }
}
