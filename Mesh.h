#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <numeric>
#include "myvectorstuff.h"
#include <iomanip>

class Mesh {
  public:

    struct Face {
      std::vector< std::vector<double> > vertices;
      std::vector<double> normal;
      double area;
      double areaCut;
      int nRayHits;
      double power;
      double Sext; //external flux incident on face
      double abs_energy_flux;
      double abs_energy_flux_previous;
      double k; //For Basko albedo scaling
      double a; //""
      double b; //""
      double previous_power;//for use in radiosity calculation
      double starting_power; //Necessary for self vf version of radiosity solver
      double Sext_previous;
      double source;
      double albedo;
      double laser_eff;
    };

    struct vfObject {
      std::string otherMesh; //Metadata to other mesh
      std::vector< std::vector<double> > F; //The view factor

        vfObject(int,int,std::string,double default_values=0.);

    };

    struct vfDatabase {
        std::vector<vfObject> content;

        vfObject operator[](std::string name){
            for(int i=0; i<content.size(); i++){
                if(name == content[i].otherMesh){return content[i];}
            }

            vfObject error(1,1,"ERROR");
            return error;
        }
    };

    std::string name;
    int numberVertices, numberFaces;
    double const_power_per_face;
    bool opaque;
    bool radiosity_converged = false; //For use in radiosity calculations
    bool self_vf;
    bool raytracing = false;
    bool basko_scaling = false;
    double min_area;
    std::vector< Face > faces;
    //Struct within Mesh to store VFs from other meshes
    vfDatabase vf_database;

    //Constructs
    Mesh (std::string, bool self_vf=false, bool opaque=false ,double const_power_per_face=0.); //prototype for constructor
    void CreateFromBlender (std::string, voxel_grid&, double albedo = 1., double basko_k = 1, double basko_a = 0, double basko_b = 0);
    //Gets / Calcs
    void viewfactor (Mesh&); //VF from this mesh to another Mesh

    //Prints
    void print_vf(std::string);
    void print_totalvf(Mesh);
    double total_area();


    void selfviewfactor();

    void print_totalselfvf();
    void block_with(Mesh&, std::vector< Mesh >);


    //For use with Voxel method (finish later)
    //std::vector<int> findinVolume(std::vector<double> origin[3], double width);
};

void saveMeshData(Mesh mesh, int timestep = 0, double realtime = 0.);
double areaFace(Mesh::Face);

#endif
