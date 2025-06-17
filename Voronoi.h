#ifndef VORONOI_HPP
#define VORONOI_HPP


#include "SaveSVG.h"

class Voronoi{
public:
    Voronoi(const std::vector<Vector> points) : points(points) {}
    void compute();
    void save(const std::string& filename) const;
    std::vector<Polygon> voronoi; 
    std::vector<Vector> points;   

protected:
    
    

    Polygon computeVoronoiCell(int idx) const;
    Polygon initializeBoundingBox() const;
    Polygon clipPolygon(const Polygon& poly, const Vector& P0, const Vector& Pi) const;


};

#endif // VORONOI_HPP