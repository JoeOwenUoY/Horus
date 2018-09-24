//
// Created by joseph on 06/05/18.
//

#ifndef HVFC_RAYTRACING_H
#define HVFC_RAYTRACING_H

#include <cmath>

#include "Mesh.h"
#include "myvectorstuff.h"

struct Ray{
    std::vector<double> r;
    std::vector<double> v;
    double power;

    Ray (std::vector<double> r, std::vector<double> v, double power){
        this->r=r;this->v=v;this->power=power; }

};

//Traces the vector of Rays to determine hits on the Mesh
//Updates the Mesh::Face.nRayHits and power (if specified).
void raytrace(std::vector<Ray>, Mesh&);

void voxel_raytrace (std::vector<Ray> rays, Mesh &mesh, voxel_grid);

void voxel_raytrace_nearest(std::vector<Ray> rays, Mesh &mesh, voxel_grid Voxel_grid);

std::vector< Ray > rays_from_txt(std::string filename, double);

#endif //HVFC_RAYTRACING_H
