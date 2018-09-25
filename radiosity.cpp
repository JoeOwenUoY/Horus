//
// Created by joseph on 07/05/18.
//

#include "radiosity.h"



//Solves the radiosity equation between two meshes
//Function will be modifying the Mesh::faces[].power variable
void solveRadiosity (Mesh &mesh1, Mesh &mesh2){

    std::string name;
    int counter = 0;
    const double tolerance = 1e-6;

    std::vector< std::vector<double> > f1, f2, fs1, fs2;

    // Including this in the radiosity solver means that the algorithm doesn't
    // have to access the very large multidimensional array from the Mesh object on
    // each array loop. Loading this in from the start produces x100 speedup!!
    f1 = mesh1.vf_database[mesh2.name].F;
    f2 = mesh2.vf_database[mesh1.name].F;
    if(mesh1.self_vf){fs1 = mesh1.vf_database[mesh1.name].F;}
    if(mesh2.self_vf){fs2 = mesh2.vf_database[mesh2.name].F;}

    //First run - to set the mesh values (allows removal of a branch in the loop (faster?))
    //Set sources
    for (int i=0; i<mesh1.faces.size(); ++i){
        mesh1.faces[i].previous_power = mesh1.faces[i].source*mesh1.faces[i].laser_eff/mesh1.faces[i].area;
    }
    for (int i=0; i<mesh2.faces.size(); ++i){
        mesh2.faces[i].previous_power = mesh2.faces[i].source*mesh2.faces[i].laser_eff/mesh2.faces[i].area;
    }

    for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

        for (unsigned int j = 0; j < mesh2.faces.size(); j++) {
            mesh1.faces[i].power += mesh1.faces[i].albedo*mesh2.faces[j].previous_power * f1[i][j];

        }
    }

    for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

        for (unsigned int j = 0; j < mesh1.faces.size(); j++) {
            mesh2.faces[i].power += mesh2.faces[i].albedo*mesh1.faces[j].previous_power * f2[i][j];

        }
    }

    //Self view factor step
    if(mesh1.self_vf) {
        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {mesh1.faces[i].previous_power = mesh1.faces[i].power;}
        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh1.faces.size(); j++) {
                mesh1.faces[i].power += mesh1.faces[i].albedo*mesh1.faces[j].previous_power * fs1[i][j];

            }
        }
    }

    //Self view factor step
    if(mesh2.self_vf){
        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {mesh2.faces[i].previous_power = mesh2.faces[i].power;}
        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh2.faces.size(); j++) {
                mesh2.faces[i].power += mesh2.faces[i].albedo*mesh2.faces[j].previous_power * fs2[i][j];

            }
        }
    }

    //Add sources
    for (int i=0; i<mesh1.faces.size(); i++){mesh1.faces[i].power += mesh1.faces[i].source*mesh1.faces[i].laser_eff/mesh1.faces[i].area;}
    for (int i=0; i<mesh2.faces.size(); i++){mesh2.faces[i].power += mesh2.faces[i].source*mesh2.faces[i].laser_eff/mesh2.faces[i].area;}

    //Main run (split from first)
    while(counter<10){

        //Set the previous mesh values
        for (unsigned int i=0; i< mesh1.faces.size(); i++){mesh1.faces[i].previous_power = mesh1.faces[i].power;
                                                           mesh1.faces[i].starting_power = mesh1.faces[i].power;
                                                           mesh1.faces[i].power = 0.;}
        for (unsigned int i=0; i< mesh2.faces.size(); i++){mesh2.faces[i].previous_power = mesh2.faces[i].power;
                                                           mesh2.faces[i].starting_power = mesh2.faces[i].power;
                                                           mesh2.faces[i].power = 0.;}

        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh2.faces.size(); j++) {

                //meshes[m].faces[j].previous_power = meshes[m].faces[j].power;
                mesh1.faces[i].power += mesh1.faces[i].albedo*mesh2.faces[j].previous_power * f1[i][j];

            }
        }


        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh1.faces.size(); j++) {

                //meshes[m].faces[j].previous_power = meshes[m].faces[j].power;
                mesh2.faces[i].power += mesh2.faces[i].albedo*mesh1.faces[j].previous_power * f2[i][j];

            }
        }

        //Self view factor step
        if(mesh1.self_vf) {
            for (unsigned int i=0; i< mesh1.faces.size(); i++){mesh1.faces[i].previous_power = mesh1.faces[i].power+
                                                               mesh1.faces[i].source*mesh1.faces[i].laser_eff/mesh1.faces[i].area;}
            for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop
                for (unsigned int j = 1; j < mesh1.faces.size(); j++) {
                    mesh1.faces[i].power += mesh1.faces[i].albedo*mesh1.faces[(i+j)%mesh1.faces.size()].previous_power * fs1[i][(i+j)%mesh1.faces.size()];

                }
            }
        }

        //Self view factor step
        if(mesh2.self_vf){
            //Move previous power setting outside so all is well before self radiosity loop
            for (unsigned int i=0; i < mesh2.faces.size();i++){mesh2.faces[i].previous_power = mesh2.faces[i].power+
                                                               mesh2.faces[i].source*mesh2.faces[i].laser_eff/mesh2.faces[i].area;}
            for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

                for (unsigned int j = 1; j < mesh2.faces.size(); j++) {
                    mesh2.faces[i].power += mesh2.faces[i].albedo*mesh2.faces[(i+j)%mesh1.faces.size()].previous_power * fs2[i][(i+j)%mesh1.faces.size()];

                }
            }
        }

        //Add sources
        for (int i=0; i<mesh1.faces.size(); i++){mesh1.faces[i].power += mesh1.faces[i].source*mesh1.faces[i].laser_eff/mesh1.faces[i].area;}
        for (int i=0; i<mesh2.faces.size(); i++){mesh2.faces[i].power += mesh2.faces[i].source*mesh2.faces[i].laser_eff/mesh2.faces[i].area;}

        mesh1.radiosity_converged = (convergence(mesh1) < tolerance)?true:false;
        //if (global_convergence(std::vector<Mesh> {mesh1, mesh2})) { return; }

        mesh2.radiosity_converged = (convergence(mesh2) < tolerance)?true:false;
        if (global_convergence(std::vector<Mesh> {mesh1, mesh2})) {counter+=1; }else{counter=0;}

    }

}




void time_dependent_solveRadiosity (Mesh &mesh1, Mesh &mesh2, std::vector< std::vector<double> >tv_record) {

    std::string name;
    int counter = 0;
    const double tolerance = 1e-6;

    std::vector<std::vector<double> > f1, f2, fs1, fs2;

    // Including this in the radiosity solver means that the algorithm doesn't
    // have to access the very large multidimensional array from the Mesh object on
    // each array loop. Loading this in from the start produces x100 speedup!!
    f1 = mesh1.vf_database[mesh2.name].F;
    f2 = mesh2.vf_database[mesh1.name].F;
    if (mesh1.self_vf) { fs1 = mesh1.vf_database[mesh1.name].F; }
    if (mesh2.self_vf) { fs2 = mesh2.vf_database[mesh2.name].F; }

    std::cout << std::scientific << "Time = "<< tv_record[0][0] << ", P_laser_in /per ray = " << tv_record[0][1] <<std::endl;

    //Inital timestep will use the arbitray initial albedo to perform the first convergence calculation.
    //Once an absorbed energy is established , then the BAsko scaling can be employed for the further timesteps below.
    //Initial timestep
    double dt = tv_record[1][0] - tv_record[0][0];

    //First run - to set the mesh values (allows removal of a branch in the loop (faster?))
    //Set sources  - now scales with time in the form of tv_record
    //Reset absorbed energy each cycle till convergence (dont want energy to runaway) not looking for energy convergence anyway
    for (int i = 0; i < mesh1.faces.size(); ++i) {

        mesh1.faces[i].previous_power =
                tv_record[0][1] * mesh1.faces[i].source * mesh1.faces[i].laser_eff * mesh1.faces[i].albedo / mesh1.faces[i].area;
        mesh1.faces[i].Sext = 0;
        mesh1.faces[i].Sext_previous = 0;
        mesh1.faces[i].abs_energy_flux = 0;

    }
    for (int i = 0; i < mesh2.faces.size(); ++i) {

        mesh2.faces[i].previous_power =
                tv_record[0][1] * mesh2.faces[i].source * mesh2.faces[i].laser_eff * mesh2.faces[i].albedo / mesh2.faces[i].area;
        mesh2.faces[i].Sext = 0;
        mesh2.faces[i].Sext_previous = 0;
        mesh2.faces[i].abs_energy_flux = 0;


    }


    for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

        for (unsigned int j = 0; j < mesh2.faces.size(); j++) {
            mesh1.faces[i].power += mesh1.faces[i].albedo * mesh2.faces[j].previous_power * f1[i][j];
        }

    }

    for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

        for (unsigned int j = 0; j < mesh1.faces.size(); j++) {
            mesh2.faces[i].power += mesh2.faces[i].albedo * mesh1.faces[j].previous_power * f2[i][j];

        }
    }

    //Self view factor step
    if (mesh1.self_vf) {
        //for (unsigned int i = 0; i < mesh1.faces.size(); i++) {mesh1.faces[i].previous_power = mesh1.faces[i].power;}
        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 1; j < mesh1.faces.size(); j++) {
                mesh1.faces[i].power += mesh1.faces[i].albedo * mesh1.faces[(i + j) % mesh1.faces.size()].previous_power * fs1[i][(i + j) % mesh1.faces.size()];

            }

        }
    }


    //Self view factor step
    if (mesh2.self_vf) {
        //for (unsigned int i = 0; i < mesh2.faces.size(); i++) {mesh2.faces[i].previous_power = mesh2.faces[i].power;}
        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh2.faces.size(); j++) {
                mesh2.faces[i].power += mesh2.faces[i].albedo * mesh2.faces[(i + j) % mesh2.faces.size()].previous_power * fs2[i][(i + j) % mesh2.faces.size()];

            }


        }
    }


    //Add sources - now scales with time in the form of tv_record
    for (int i = 0; i < mesh1.faces.size(); i++) { mesh1.faces[i].power += tv_record[0][1] * mesh1.faces[i].source* mesh1.faces[i].albedo *
                                                                           mesh1.faces[i].laser_eff / mesh1.faces[i].area;
    }
    for (int i = 0; i < mesh2.faces.size(); i++) { mesh2.faces[i].power += tv_record[0][1] * mesh2.faces[i].source * mesh2.faces[i].albedo *
                                                                           mesh2.faces[i].laser_eff / mesh2.faces[i].area;
    }


    //Main run (split from first)
    while (counter < 5) {

        //Set the previous mesh values
        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {
            mesh1.faces[i].previous_power = mesh1.faces[i].power;
            mesh1.faces[i].starting_power = mesh1.faces[i].power;
            mesh1.faces[i].power = 0.;
            mesh1.faces[i].abs_energy_flux = 0;
            mesh1.faces[i].Sext_previous = mesh1.faces[i].Sext;
            mesh1.faces[i].Sext = 0;
        }
        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {
            mesh2.faces[i].previous_power = mesh2.faces[i].power;
            mesh2.faces[i].starting_power = mesh2.faces[i].power;
            mesh2.faces[i].power = 0.;
            mesh2.faces[i].abs_energy_flux = 0;
            mesh2.faces[i].Sext_previous = mesh2.faces[i].Sext;
            mesh2.faces[i].Sext = 0;
        }

        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh2.faces.size(); j++) {

                mesh1.faces[i].Sext += mesh2.faces[j].previous_power * f1[i][j];


            }
        }


        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh1.faces.size(); j++) {

                //meshes[m].faces[j].previous_power = meshes[m].faces[j].power;
                mesh2.faces[i].Sext += mesh1.faces[j].previous_power * f2[i][j];

            }
        }

        //Self view factor step
        if (mesh1.self_vf) {
            //                for (unsigned int i = 0; i < mesh1.faces.size(); i++) {
            //                    mesh1.faces[i].previous_power = mesh1.faces[i].power +
            //                                                    mesh1.faces[i].source * mesh1.faces[i].laser_eff /
            //                                                    mesh1.faces[i].area;
            //                }
            for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop
                for (unsigned int j = 1; j < mesh1.faces.size(); j++) {
                    mesh1.faces[i].Sext += mesh1.faces[(i + j) % mesh1.faces.size()].previous_power *
                                           fs1[i][(i + j) % mesh1.faces.size()];

                }
            }
        }

        //Self view factor step
        if (mesh2.self_vf) {
            //                //Move previous power setting outside so all is well before self radiosity loop
            //                for (unsigned int i = 0; i < mesh2.faces.size(); i++) {
            //                    mesh2.faces[i].previous_power = mesh2.faces[i].power +
            //                                                    mesh2.faces[i].source * mesh2.faces[i].laser_eff /
            //                                                    mesh2.faces[i].area;
            //                }
            for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

                for (unsigned int j = 1; j < mesh2.faces.size(); j++) {
                    mesh2.faces[i].Sext += mesh2.faces[(i + j) % mesh2.faces.size()].previous_power *
                                           fs2[i][(i + j) % mesh2.faces.size()];


                }

            }
        }


        mesh1.radiosity_converged = (convergence_Sext(mesh1) < tolerance) ? true : false;
        //if (global_convergence(std::vector<Mesh> {mesh1, mesh2})) { return; }
        //std::cout << "Conv 1 " <<convergence(mesh1) << "tol. "<<tolerance<<std::endl;
        mesh2.radiosity_converged = (convergence_Sext(mesh2) < tolerance) ? true : false;
        if (global_convergence(std::vector<Mesh>{mesh1, mesh2})) { counter += 1; } else { counter = 0; }
        //std::cout << "Conv 2. " <<convergence(mesh2) << "tol. "<<tolerance<<std::endl;

    }

    //Calc the reflected radiation flux - Basko scaling
    //Mesh 1
    for (int i = 0; i < mesh1.faces.size(); i++) {

        mesh1.faces[i].Sext += tv_record[0][1] * mesh1.faces[i].source * mesh1.faces[i].laser_eff / mesh1.faces[i].area;

        mesh1.faces[i].power = mesh1.faces[i].albedo * mesh1.faces[i].Sext;
        //Calc absorbed energy
        mesh1.faces[i].abs_energy_flux = mesh1.faces[i].Sext * (1.0 - mesh1.faces[i].albedo) * dt;

        mesh1.faces[i].albedo = flux_convergence(mesh1.faces[i].Sext*(1e-14)/(1e4), mesh1.faces[i].k,
                                                 mesh1.faces[i].abs_energy_flux*(1e-6)/(1e4),
                                                 mesh1.faces[i].a, mesh1.faces[i].b, mesh1.faces[i].albedo);




    }
    //Mesh 2
    for (int i = 0; i < mesh2.faces.size(); i++) {

        mesh2.faces[i].Sext += tv_record[0][1] * mesh2.faces[i].source * mesh2.faces[i].laser_eff / mesh2.faces[i].area;

        mesh2.faces[i].power = mesh2.faces[i].albedo * mesh2.faces[i].Sext;
        //Calc absorbed energy
        mesh2.faces[i].abs_energy_flux += dt*(1.0 - mesh2.faces[i].albedo) * (mesh2.faces[i].Sext);

        mesh2.faces[i].albedo = flux_convergence(mesh2.faces[i].Sext*(1e-14)/(1e4), mesh2.faces[i].k,
                                                 mesh2.faces[i].abs_energy_flux*(1e-6)/(1e4),
                                                 mesh2.faces[i].a, mesh2.faces[i].b, mesh2.faces[i].albedo);

    }



    saveMeshData(mesh1,1);
    saveMeshData(mesh2,1);

    // ########### START OF TIMESTEP ###############

    for (int timestep = 1; timestep < tv_record.size(); timestep++) {
        std::cout << std::scientific << "Time = "<< tv_record[timestep][0] << ", P_laser_in / per ray = " << tv_record[timestep][1] <<std::endl;

        counter = 0;

        dt = (timestep != tv_record.size() - 1) ? tv_record[timestep + 1][0] - tv_record[timestep][0] : 0.;

        //First run - to set the mesh values (allows removal of a branch in the loop (faster?))
        //Set sources  - now scales with time in the form of tv_record
        for (int i = 0; i < mesh1.faces.size(); ++i) {
            mesh1.faces[i].previous_power =
                    tv_record[timestep][1] * mesh1.faces[i].source * mesh1.faces[i].laser_eff * mesh1.faces[i].albedo/mesh1.faces[i].area;
            mesh1.faces[i].Sext = tv_record[timestep][1] * mesh1.faces[i].source * mesh1.faces[i].laser_eff /mesh1.faces[i].area;
            mesh1.faces[i].Sext_previous = 0;
            mesh1.faces[i].abs_energy_flux_previous = mesh1.faces[i].abs_energy_flux;


        }
        for (int i = 0; i < mesh2.faces.size(); ++i) {
            mesh2.faces[i].previous_power =
                    tv_record[timestep][1] * mesh2.faces[i].source * mesh2.faces[i].laser_eff * mesh2.faces[i].albedo/mesh2.faces[i].area;
            mesh2.faces[i].Sext = tv_record[timestep][1] * mesh2.faces[i].source * mesh2.faces[i].laser_eff /mesh2.faces[i].area;
            mesh2.faces[i].Sext_previous = 0;
            mesh2.faces[i].abs_energy_flux_previous = mesh2.faces[i].abs_energy_flux;

        }

        for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh2.faces.size(); j++) {
                mesh1.faces[i].Sext += mesh2.faces[j].previous_power * f1[i][j];
            }
        }

        for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

            for (unsigned int j = 0; j < mesh1.faces.size(); j++) {
                mesh2.faces[i].Sext += mesh1.faces[j].previous_power * f2[i][j];

            }
        }

        //Self view factor step
        if (mesh1.self_vf) {
            //                for (unsigned int i = 0; i < mesh1.faces.size(); i++) {
            //                    mesh1.faces[i].previous_power = mesh1.faces[i].power +
            //                                                    mesh1.faces[i].source * mesh1.faces[i].laser_eff /
            //                                                    mesh1.faces[i].area;
            //                }
            for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop
                for (unsigned int j = 1; j < mesh1.faces.size(); j++) {
                    mesh1.faces[i].Sext += mesh1.faces[(i + j) % mesh1.faces.size()].previous_power *
                                           fs1[i][(i + j) % mesh1.faces.size()];

                }
            }
        }

        //Self view factor step
        if (mesh2.self_vf) {
            //                //Move previous power setting outside so all is well before self radiosity loop
            //                for (unsigned int i = 0; i < mesh2.faces.size(); i++) {
            //                    mesh2.faces[i].previous_power = mesh2.faces[i].power +
            //                                                    mesh2.faces[i].source * mesh2.faces[i].laser_eff /
            //                                                    mesh2.faces[i].area;
            //                }
            for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

                for (unsigned int j = 1; j < mesh2.faces.size(); j++) {
                    mesh2.faces[i].Sext += mesh2.faces[(i + j) % mesh2.faces.size()].previous_power *
                                           fs2[i][(i + j) % mesh2.faces.size()];

                }
            }
        }


        //Calc the reflected radiation flux - Basko scaling
        //Mesh 1
        for (int i = 0; i < mesh1.faces.size(); i++) {

//            //calc abledo based on external flux incident (not the lasers though, these are accounted for by
//            // the linear scaling conversion parameter (correct?))) S is in 10^14 W/cc, E_a is in MJ/cc
//            mesh1.faces[i].albedo = flux_convergence(mesh1.faces[i].Sext*(1e-14)/(mesh1.faces[i].area*1e4), mesh1.faces[i].k,
//                                                     mesh1.faces[i].abs_energy_flux*(1e-6)/(mesh1.faces[i].area*1e4),
//                                                     mesh1.faces[i].a, mesh1.faces[i].b, mesh1.faces[i].albedo);
            mesh1.faces[i].power = mesh1.faces[i].albedo * mesh1.faces[i].Sext;
//            mesh1.faces[i].power += tv_record[timestep][1] * mesh1.faces[i].source *
//                                    mesh1.faces[i].laser_eff; //Add sources
//            //Calc absorbed energy
//            mesh1.faces[i].abs_energy_flux += dt*((1.0 - mesh1.faces[i].albedo) * mesh1.faces[i].Sext);
        }
        //Mesh 2
        for (int i = 0; i < mesh2.faces.size(); i++) {

//            mesh2.faces[i].albedo = flux_convergence(mesh2.faces[i].Sext*(1e-14)/(mesh2.faces[i].area*1e4), mesh2.faces[i].k,
//                                                     mesh2.faces[i].abs_energy_flux*(1e-6)/(mesh2.faces[i].area*1e4),
//                                                     mesh2.faces[i].a, mesh2.faces[i].b, mesh2.faces[i].albedo);
            mesh2.faces[i].power = mesh2.faces[i].albedo * mesh2.faces[i].Sext;
//            mesh2.faces[i].power += tv_record[timestep][1] * mesh2.faces[i].source *
//                                    mesh2.faces[i].laser_eff; //Add sources
//            //Calc absorbed energy
//            mesh2.faces[i].abs_energy_flux += dt*((1.0 - mesh2.faces[i].albedo) * mesh2.faces[i].Sext);
        }




        //Main run (split from first)
        while (counter < 5) {

            //Set the previous mesh values
            for (unsigned int i = 0; i < mesh1.faces.size(); i++) {
                mesh1.faces[i].previous_power = mesh1.faces[i].power;
                mesh1.faces[i].starting_power = mesh1.faces[i].power;

               // mesh1.faces[i].abs_energy_flux = mesh1.faces[i].abs_energy_flux_previous;
                mesh1.faces[i].Sext_previous = mesh1.faces[i].Sext;
                mesh1.faces[i].Sext = tv_record[timestep][1] * mesh1.faces[i].source * mesh1.faces[i].laser_eff/mesh1.faces[i].area;;
            }
            for (unsigned int i = 0; i < mesh2.faces.size(); i++) {
                mesh2.faces[i].previous_power = mesh2.faces[i].power;
                mesh2.faces[i].starting_power = mesh2.faces[i].power;
                mesh2.faces[i].Sext_previous = mesh2.faces[i].Sext;

               // mesh2.faces[i].abs_energy_flux = mesh2.faces[i].abs_energy_flux_previous;
                mesh2.faces[i].Sext = tv_record[timestep][1] * mesh2.faces[i].source * mesh2.faces[i].laser_eff/mesh2.faces[i].area;;
            }


            

            for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop

                for (unsigned int j = 0; j < mesh2.faces.size(); j++) {

                    mesh1.faces[i].Sext += mesh2.faces[j].previous_power * f1[i][j];

                }
            }


            for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

                for (unsigned int j = 0; j < mesh1.faces.size(); j++) {

                    //meshes[m].faces[j].previous_power = meshes[m].faces[j].power;
                    mesh2.faces[i].Sext += mesh1.faces[j].previous_power * f2[i][j];

                }
            }

            //Self view factor step
            if (mesh1.self_vf) {
                //                for (unsigned int i = 0; i < mesh1.faces.size(); i++) {
                //                    mesh1.faces[i].previous_power = mesh1.faces[i].power +
                //                                                    mesh1.faces[i].source * mesh1.faces[i].laser_eff /
                //                                                    mesh1.faces[i].area;
                //                }
                for (unsigned int i = 0; i < mesh1.faces.size(); i++) {//first mesh faces loop
                    for (unsigned int j = 1; j < mesh1.faces.size(); j++) {
                        mesh1.faces[i].Sext += mesh1.faces[(i + j) % mesh1.faces.size()].previous_power *
                                               fs1[i][(i + j) % mesh1.faces.size()];

                    }
                }
            }

            //Self view factor step
            if (mesh2.self_vf) {
                //                //Move previous power setting outside so all is well before self radiosity loop
                //                for (unsigned int i = 0; i < mesh2.faces.size(); i++) {
                //                    mesh2.faces[i].previous_power = mesh2.faces[i].power +
                //                                                    mesh2.faces[i].source * mesh2.faces[i].laser_eff /
                //                                                    mesh2.faces[i].area;
                //                }
                for (unsigned int i = 0; i < mesh2.faces.size(); i++) {//first mesh faces loop

                    for (unsigned int j = 1; j < mesh2.faces.size(); j++) {
                        mesh2.faces[i].Sext += mesh2.faces[(i + j) % mesh2.faces.size()].previous_power *
                                               fs2[i][(i + j) % mesh2.faces.size()];

                    }
                }
            }




            mesh1.radiosity_converged = (convergence_Sext(mesh1) < tolerance) ? true : false;
            //if (global_convergence(std::vector<Mesh> {mesh1, mesh2})) { return; }
            //std::cout << "Conv 1 " <<convergence(mesh1) << "tol. "<<tolerance<<std::endl;

            mesh2.radiosity_converged = (convergence_Sext(mesh2) < tolerance) ? true : false;
            if (global_convergence(std::vector<Mesh>{mesh1, mesh2})) { counter += 1; } else { counter = 0; }
            //std::cout << "Conv 2 " <<convergence(mesh2) << "tol. "<<tolerance<<std::endl;
        }

        //Calc the reflected radiation flux - Basko scaling
        //Mesh 1
        for (int i = 0; i < mesh1.faces.size(); i++) {

            mesh1.faces[i].power = mesh1.faces[i].albedo * mesh1.faces[i].Sext;
            //Calc absorbed energy
            mesh1.faces[i].abs_energy_flux += mesh1.faces[i].Sext * (1.0 - mesh1.faces[i].albedo) * dt;

           mesh1.faces[i].albedo = flux_convergence(mesh1.faces[i].Sext*(1e-14)/(1e4), mesh1.faces[i].k,
                                                  mesh1.faces[i].abs_energy_flux*(1e-6)/(1e4),
                                                     mesh1.faces[i].a, mesh1.faces[i].b, mesh1.faces[i].albedo);

        }
        //Mesh 2
        for (int i = 0; i < mesh2.faces.size(); i++) {

            mesh2.faces[i].power = mesh2.faces[i].albedo * mesh2.faces[i].Sext;
            //Calc absorbed energy
            mesh2.faces[i].abs_energy_flux += mesh2.faces[i].Sext * (1.0 - mesh2.faces[i].albedo) * dt;


            mesh2.faces[i].albedo = flux_convergence(mesh2.faces[i].Sext*(1e-14)/(1e4), mesh2.faces[i].k,
                                                     mesh2.faces[i].abs_energy_flux*(1e-6)/(1e4),
                                                     mesh2.faces[i].a, mesh2.faces[i].b, mesh2.faces[i].albedo);
        }


        saveMeshData(mesh1, timestep+1);
        saveMeshData(mesh2, timestep+1);

        // ########### END OF TIMESTEP ###############


    }

}