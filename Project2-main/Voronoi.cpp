#include "Voronoi.h"


void Voronoi::compute() {
    voronoi.resize(points.size());
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        voronoi[i] = computeVoronoiCell(i);
    }
}

void Voronoi::save(const std::string& filename) const {
    save_svg(voronoi, filename, "white");
}

Polygon Voronoi::computeVoronoiCell(int idx) const {
    Polygon cell;
    cell.vertices = {
        Vector(0, 0, 0),
        Vector(0, 1, 0),
        Vector(1, 1, 0),
        Vector(1, 0, 0)
    };

    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        if (i != idx){
            cell = clipPolygon(cell, points[idx], points[i]);
        }
    }

    return cell;
}

Polygon Voronoi::clipPolygon(const Polygon& poly, const Vector& P0, const Vector& Pi) const {
    Polygon result;
    Vector M = (P0 + Pi) / 2;

    for (size_t i = 0; i < poly.vertices.size(); ++i) {
        Vector A = (i == 0) ? poly.vertices.back() : poly.vertices[i - 1];
        const Vector& B = poly.vertices[i];

        double t = dot(M - A, Pi - P0) / dot(B - A, Pi - P0);
        Vector P = A + t * (B - A);

        double norm2B_P0 = (B - P0).norm2();
        double norm2B_Pi = (B - Pi).norm2();
        double norm2A_P0 = (A - P0).norm2();
        double norm2A_Pi = (A - Pi).norm2();

        if (norm2B_P0 < norm2B_Pi) {
            if (norm2A_P0 > norm2A_Pi) {
                result.vertices.push_back(P);
            }
            result.vertices.push_back(B);
        } else if (norm2A_P0 < norm2A_Pi) {
            result.vertices.push_back(P);
        }
    }

    return result;
}
