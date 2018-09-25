//
// Created by joseph on 07/05/18.
//

#ifndef HVFC_RADIOSITY_H
#define HVFC_RADIOSITY_H

#include "myvectorstuff.h"
#include "Mesh.h"
#include "time_dependence.h"

inline double convergence(Mesh mesh) {
    double max_difference = 0;
    for (int i = 0; i<mesh.faces.size(); i++){
        max_difference = (max_difference < fabs(mesh.faces[i].starting_power-mesh.faces[i].power))?fabs(mesh.faces[i].starting_power-mesh.faces[i].power):
                         max_difference;
    }
    //std::cout << mesh.name<<": "<<max_difference <<std::endl;
    return max_difference;
}

inline double convergence_Sext(Mesh mesh) {
    double max_difference = 0;
    for (int i = 0; i<mesh.faces.size(); i++){
        max_difference = (max_difference < fabs(mesh.faces[i].Sext_previous-mesh.faces[i].Sext))?fabs(mesh.faces[i].Sext_previous-mesh.faces[i].Sext):
                         max_difference;
    }
    //std::cout << mesh.name<<": "<<max_difference <<std::endl;
    return max_difference;
}


inline bool global_convergence (std::vector<Mesh> scene) {
    for (int i=0; i<scene.size(); i++){
        if(!scene[i].radiosity_converged){return false;}
    }
    return true;
}

void solveRadiosity (Mesh&, Mesh&);
void time_dependent_solveRadiosity (Mesh&, Mesh&, std::vector< std::vector<double> >);
inline bool global_convergence(std::vector< Mesh >);
inline double convergence (std::vector<double>, Mesh);

#endif //HVFC_RADIOSITY_H
