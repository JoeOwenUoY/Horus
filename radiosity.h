//
// Created by joseph on 07/05/18.
//

#ifndef HVFC_RADIOSITY_H
#define HVFC_RADIOSITY_H

#include "myvectorstuff.h"
#include "headers/Mesh.h"
#include "headers/time_dependence.h"

void solveRadiosity (Mesh&, Mesh&);
void time_dependent_solveRadiosity (Mesh&, Mesh&, std::vector< std::vector<double> >);
inline bool global_convergence(std::vector< Mesh >);
inline double convergence (std::vector<double>, Mesh);

#endif //HVFC_RADIOSITY_H
