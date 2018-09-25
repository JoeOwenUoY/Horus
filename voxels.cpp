//
// Created by jmo510 on 06/06/18.
//

#include "voxels.h"

voxel_grid create_voxel_grid(std::vector<double> grid_center, double x_width, double y_width, double z_width, int nVoxels) {

    voxel_grid Voxel_grid;
    Voxel_grid.x_width=x_width;
    Voxel_grid.y_width=y_width;
    Voxel_grid.z_width=z_width;
    Voxel_grid.center=grid_center;

    double x_density = double(nVoxels)/x_width;
    double y_density = double(nVoxels)/x_width;
    double z_density = double(nVoxels)/x_width;

    Voxel_grid.x_cell_size = 1.0/x_density;
    Voxel_grid.y_cell_size = 1.0/y_density;
    Voxel_grid.z_cell_size = 1.0/z_density;

    std::cout << x_density<< " "<< y_density<< " "<<z_density<<std::endl;

    std::vector<double> starting_center = grid_center - std::vector<double>{x_width/2.0,y_width/2.0,z_width/2.0};

    for (int i=0; i<int(x_density*x_width); i++) {
        for (int j = 0; j < int(y_density*y_width); j++) {
            for (int k = 0; k < int(z_density*z_width); k++) {


                Voxel_grid.grid.push_back(voxel(starting_center + std::vector<double>{double(i)/x_density, double(j)/y_density, double(k)/z_density},
                                           1.0/x_density, 1.0/y_density, 1.0/z_density ));

            }

        }
    }
    std::cout << "# voxels = "<< Voxel_grid.grid.size()<<std::endl;

    return Voxel_grid;
}



void save_voxel_grid(voxel_grid Voxel_grid) {
    std::ofstream save ("./voxelgrid.save");

    for (auto voxel : Voxel_grid.grid) {
        if (voxel.content.size() > 0) {
            save << voxel.center[0] << " " << voxel.center[1] << " " << voxel.center[2] << std::endl;
        }
    }
    return;
}
