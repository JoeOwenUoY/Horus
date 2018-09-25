//
// Created by joseph on 06/05/18.
//

#include <limits>
#include "raytracing.h"




////For the string input
//template <typename MyType>
//std::vector<MyType> String_to_Vector (std::string line){
//    std::istringstream iss(line);
//    std::vector<MyType> numbers;
//    MyType T;
//    while (iss >> T ){
//        numbers.push_back( T );
//    }
//    return numbers;
//}


void raytrace (std::vector<Ray> rays, Mesh &mesh){

    double t;
    std::vector<double> r (3,0.);
    std::vector<double> p (3,0.);
    std::vector<double> d (3,0.);//vector to solution
    double alpha, beta;
    const int nSegments = 5;
    double tolerance;
    int hitcounter = 0;


    for(unsigned int i = 0; i<rays.size(); i++){
        for(unsigned int j = 0; j<mesh.faces.size(); j++){
            tolerance = sqrt(mesh.faces[j].area)/(nSegments);

            if(dot(rays[i].v,mesh.faces[j].normal)<0){
                t = dot((mesh.faces[j].vertices[0]+
                        mesh.faces[j].vertices[1]+
                        mesh.faces[j].vertices[2])/3.0 - rays[i].r,
                        mesh.faces[j].normal)/dot(rays[i].v,mesh.faces[j].normal);

                p = t*rays[i].v + rays[i].r;

                //Check if the point lies within the face
                for(int a = 0; a<=nSegments; a++) {
                    alpha = 1.0 * a / (1.0 * nSegments);
                    for(int b = 0; b<=nSegments-a; b++) {
                        beta = 1.0 * b / (1.0 * nSegments);
                        d = (alpha * (mesh.faces[j].vertices[1] - mesh.faces[j].vertices[0]) +
                            beta * (mesh.faces[j].vertices[2] - mesh.faces[j].vertices[0]))
                            -(p-mesh.faces[j].vertices[0]);

                        if (magnitude(d) < tolerance) {

                            mesh.faces[j].nRayHits += 1;
                            mesh.faces[j].source += rays[i].power;
                            hitcounter += 1;
                            goto ray_hit;
                        }
                    }
                }

                for(int b = 0; b<=nSegments; b++) {
                    beta = 1.0 * b / (1.0 * nSegments);
                    for(int a = 0; a<=nSegments-a; a++) {
                        alpha = 1.0 * a / (1.0 * nSegments);
                        d = (alpha * (mesh.faces[j].vertices[1] - mesh.faces[j].vertices[0]) +
                             beta * (mesh.faces[j].vertices[2] - mesh.faces[j].vertices[0]))
                            -(p-mesh.faces[j].vertices[0]);

                        if (magnitude(d) < tolerance) {

                            mesh.faces[j].nRayHits += 1;
                            mesh.faces[j].source += rays[i].power;
                            hitcounter += 1;
                            goto ray_hit;
                        }
                    }
                }




                    }
                }
        ray_hit:

        continue;//GOTO for when a ray hits
            }

            std::cout<<"# Hits = "<<hitcounter<<std::endl;

}

std::vector<Ray> rays_from_txt(std::string filename, double total_power) {
    std::ifstream myfile (filename);
    std::string line;
    std::vector<Ray> rays;


    if (myfile.is_open()){
        while( getline(myfile, line) ){
            rays.push_back(Ray (std::vector<double> {String_to_Vector<double>(line)[0],String_to_Vector<double>(line)[1],String_to_Vector<double>(line)[2]},
            std::vector<double> {String_to_Vector<double>(line)[3],
                    String_to_Vector<double>(line)[4],String_to_Vector<double>(line)[5]}, String_to_Vector<double>(line)[6]));
        }
    }else {std::cout<<"Could not open ray file!"<<std::endl;}
    myfile.close();

    for (int i=0; i<rays.size(); i++){
        rays[i].power=1;
    }



    return rays;
}



void voxel_raytrace(std::vector<Ray> rays, Mesh &mesh, voxel_grid Voxel_grid) {

    double t;
    std::vector<double> r (3,0.);
    std::vector<double> p (3,0.);
    std::vector<double> d (3,0.);//vector to solution
    double alpha, beta;
    const int nSegments = 15;
    double tolerance;
    int hitcounter = 0;

    //The ray actually propagates in this method
    std::vector<double> ray_location (3,0.);
    int voxel_index; //Will be index of nearest voxel
    bool ray_out_of_bounds;

    double max_ray_length = sqrt(Voxel_grid.x_width*Voxel_grid.x_width +
                             Voxel_grid.y_width*Voxel_grid.y_width +
                             Voxel_grid.z_width*Voxel_grid.z_width);

    int step_density = 5;
    int max_tsteps = step_density*(pow(Voxel_grid.grid.size(),1.0/3.0));

    for(unsigned int i = 0; i<rays.size(); i++){
        //Loop over the rays
        for (int timestep=0; timestep<max_tsteps; timestep++){
            //Propagate the ray
            ray_location = (max_ray_length*timestep/max_tsteps)*rays[i].v + rays[i].r;
            voxel_index=nearest_point2point(ray_location, Voxel_grid.grid);

            //Loop over faces contained in the voxel instead
            for(voxel::member Member : Voxel_grid.grid[voxel_index].content){
                //If the mesh have a face in the voxel
                if(mesh.name==Member.mesh_name){
                    //The rest is the usual raytracing algorithm as in the normal routine
                    //but with only faces in the voxel looped over
                    tolerance = sqrt(mesh.faces[Member.face_index].area)/(nSegments);

                    if(dot(rays[i].v,mesh.faces[Member.face_index].normal)<0){
                        t = dot((mesh.faces[Member.face_index].vertices[0]+
                                 mesh.faces[Member.face_index].vertices[1]+
                                 mesh.faces[Member.face_index].vertices[2])/3.0 - rays[i].r,
                                mesh.faces[Member.face_index].normal)/dot(rays[i].v,mesh.faces[Member.face_index].normal);

                        p = t*rays[i].v + rays[i].r;

                        //Check if the point lies within the face
                        for(int a = 0; a<=nSegments; a++) {
                            alpha = 1.0 * a / (1.0 * nSegments);
                            for(int b = 0; b<=nSegments-a; b++) {
                                beta = 1.0 * b / (1.0 * nSegments);
                                d = (alpha * (mesh.faces[Member.face_index].vertices[1] - mesh.faces[Member.face_index].vertices[0]) +
                                     beta * (mesh.faces[Member.face_index].vertices[2] - mesh.faces[Member.face_index].vertices[0]))
                                    -(p-mesh.faces[Member.face_index].vertices[0]);

                                if (magnitude(d) < tolerance) {

                                    mesh.faces[Member.face_index].nRayHits += 1;
                                    mesh.faces[Member.face_index].source += rays[i].power;
                                    hitcounter += 1;
                                    goto ray_hit;
                                }
                            }
                        }

                        for(int b = 0; b<=nSegments; b++) {
                            beta = 1.0 * b / (1.0 * nSegments);
                            for(int a = 0; a<=nSegments-a; a++) {
                                alpha = 1.0 * a / (1.0 * nSegments);
                                d = (alpha * (mesh.faces[Member.face_index].vertices[1] - mesh.faces[Member.face_index].vertices[0]) +
                                     beta * (mesh.faces[Member.face_index].vertices[2] - mesh.faces[Member.face_index].vertices[0]))
                                    -(p-mesh.faces[Member.face_index].vertices[0]);

                                if (magnitude(d) < tolerance) {

                                    mesh.faces[Member.face_index].nRayHits += 1;
                                    mesh.faces[Member.face_index].source += rays[i].power;
                                    hitcounter += 1;
                                    goto ray_hit;
                                }
                            }
                        }




                    }
                }
            }
        }

        ray_hit:

        continue;//GOTO for when a ray hits
    }

    std::cout<<"# Hits = "<<hitcounter<<std::endl;
}

// Uses a faster raytrace hit detection, perhaps less accurate
//A version of voxell raytrace but with a (hopefully) better chance of deciding which face was the one hit.
//Looks for the face with the closed approach to the ray that would allow for contact, defined by the initial
//value of 'nearest', which is the same as 'tolerance' in the origin fuction
void voxel_raytrace_nearest(std::vector<Ray> rays, Mesh &mesh, voxel_grid Voxel_grid) {

    double t;
    std::vector<double> r (3,0.);
    std::vector<double> p (3,0.);
    std::vector<double> d (3,0.);//vector to solution
    double alpha, beta;
    const int nSegments = 15;
    int hitcounter = 0;
    int nearestid;
    double nearest;
    std::string nearest_name;
    //The ray actually propagates in this method
    std::vector<double> ray_location (3,0.);
    int voxel_index; //Will be index of nearest voxel
    bool ray_out_of_bounds;

    double max_ray_length = sqrt(Voxel_grid.x_width*Voxel_grid.x_width +
                                 Voxel_grid.y_width*Voxel_grid.y_width +
                                 Voxel_grid.z_width*Voxel_grid.z_width);

    std::vector<double> centroid (3,0.);

    //Important parameter MAYBE PUT SOMEWHERE ELSE LATER

    int step_density = 5;
    int max_tsteps = step_density*(pow(Voxel_grid.grid.size(),1.0/3.0));

    for(unsigned int i = 0; i<rays.size(); i++) {

        //Reset hit parameters
        nearestid = -1;
        nearest = pow(mesh.min_area,0.5); //Effective distance of closest approach for a ray


        //Loop over the rays
        for (int timestep = 0; timestep < max_tsteps; timestep++) {
            //Propagate the ray
            ray_location = (max_ray_length * timestep / max_tsteps) * rays[i].v + rays[i].r;
            voxel_index = nearest_point2point(ray_location, Voxel_grid.grid);

            //Loop over faces contained in the voxel instead
            for (voxel::member Member : Voxel_grid.grid[voxel_index].content) {
                //If the mesh have a face in the voxel
                if (mesh.name == Member.mesh_name) {
                    //The rest is the usual raytracing algorithm as in the normal routine
                    //but with only faces in the voxel looped over


                    if (dot(rays[i].v, mesh.faces[Member.face_index].normal) < 0) {
                        t = dot((mesh.faces[Member.face_index].vertices[0] +
                                 mesh.faces[Member.face_index].vertices[1] +
                                 mesh.faces[Member.face_index].vertices[2]) / 3.0 - rays[i].r,
                                mesh.faces[Member.face_index].normal) /
                            dot(rays[i].v, mesh.faces[Member.face_index].normal);

                        p = t * rays[i].v + rays[i].r;
                        centroid = (mesh.faces[Member.face_index].vertices[0] +
                                    mesh.faces[Member.face_index].vertices[1] +
                                    mesh.faces[Member.face_index].vertices[2]) / 3.0;
                        //Check if it is some distance from the face centroid

                        if (magnitude(p - centroid) < nearest) {

                            nearest = magnitude(p - centroid);
                            nearestid = Member.face_index;

                        }

                    }


                }
            }
        }


        if (nearestid != -1) {
            mesh.faces[nearestid].nRayHits += 1;
            mesh.faces[nearestid].source += rays[i].power;
            hitcounter += 1;
        }

    }

    std::cout<<"# Hits = "<<hitcounter<<std::endl;
}


