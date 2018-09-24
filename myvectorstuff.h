#ifndef MYVECTORSTUFF_H
#define MYVECTORSTUFF_H

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>

#include "voxels.h"

double magnitude(std::vector<double> x);
std::vector<double> midpoint(std::vector<double>, std::vector<double>,
                             std::vector<double>, std::vector<double>);

int face_visible(std::vector<double>,std::vector<double>);
double dot(std::vector<double>,std::vector<double>);

template <typename T>
void transpose(std::vector< std::vector<T> >);

void rotate(std::vector<double> &vector, std::vector<double> axis, double angle );


int nearest_point2point (std::vector<double> point, std::vector<voxel> voxels);
#endif
