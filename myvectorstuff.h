#ifndef MYVECTORSTUFF_H
#define MYVECTORSTUFF_H

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include "voxels.h"

class voxel;

//Overload + operator for vectors
template <typename T>
inline std::vector<T> operator+(std::vector<T> x,std::vector<T> y){
    if(x.size()==y.size()){
        std::vector<T> z (x.size(),0.);
        for(int i=0; i<x.size(); i++){
            z[i] = x[i] + y[i];
        }
        return z;
    }else{
        std::cout << "Vectors x and y must be the same size."<<std::endl;
        std::exit(EXIT_FAILURE);
    }
}

template <typename T>
inline std::vector<T> operator-(std::vector<T> x,std::vector<T> y){
    if(x.size()==y.size()){
        std::vector<T> z (x.size(),0.);
        for(int i=0; i<x.size(); i++){
            z[i] = x[i] - y[i];
        }
        return z;
    }else{
        std::cout << "Vectors x and y must be the same size."<<std::endl;
        std::exit(EXIT_FAILURE);
    }
}

template <typename T>
inline std::vector<T> operator/(std::vector<T> x, double a){
    std::vector<T> z (x.size(),0.);
    for(int i=0; i<x.size(); i++){
        z[i] = x[i]/a;
    }
    return z;
}

template <typename T>
inline std::vector<T> operator*( double a, std::vector<T> x){
    std::vector<T> z (x.size(),0.);
    for(int i=0; i<x.size(); i++){
        z[i] = x[i]*a;
    }
    return z;
}


//matrix mul a vector//
template <typename T>
inline std::vector<T> operator*(std::vector< std::vector<T> > M, std::vector<T> v){
    std::vector<T> result (v.size(),0.);
    for(int j=0; j<v.size(); j++){
        for (int i=0; i<v.size(); i++){
            result[j] += M[j][i]*v[i];
        }
    }
    return result;
}

//For the string input
template <typename MyType>
std::vector<MyType> String_to_Vector (std::string line){
    std::istringstream iss(line);
    std::vector<MyType> numbers;
    MyType T;
    while (iss >> T ){
        numbers.push_back( T );
    }
    return numbers;
}

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
