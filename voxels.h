//
// Created by jmo510 on 06/06/18.
//

#ifndef HVFC_VOXELS_H
#define HVFC_VOXELS_H

#include <fstream>
#include "myvectorstuff.h"


//Quantize 3D space into 'voxels'
//A single voxel has its own dimensions and a list of meshes who have faces inside it as well
//as the index of the faces for each mesh.

class voxel {

public:

    struct member{
        std::string mesh_name;
        int face_index;
    };

    std::vector<double> center;
    double x_width;
    double y_width;
    double z_width;
    std::vector<member> content;

    voxel (std::vector<double> center, double x_width, double y_width, double z_width){
        this->center=center;
        this->x_width = x_width;
        this->y_width = y_width;
        this->z_width = z_width;
    }

};

struct voxel_grid {
    //The actual grid of voxels
    std::vector<voxel> grid;
    //The properties of the grid
    double x_width;
    double y_width;
    double z_width;
    std::vector<double> center;
    //Component voxel properties
    double x_cell_size;
    double y_cell_size;
    double z_cell_size;

};


voxel_grid create_voxel_grid (std::vector<double> grid_center, double x_width, double y_width, double z_width, int nVoxels);
void save_voxel_grid (voxel_grid Voxel_grid);
#endif //HVFC_VOXELS_H
