#include "myvectorstuff.h"



double dot(std::vector<double> x, std::vector<double> y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

double magnitude(std::vector<double> x) {
  return std::sqrt(dot(x,x));
}

std::vector<double> midpoint(std::vector<double> x1, std::vector<double> x2,
                             std::vector<double> y1, std::vector<double> y2){
  return std::vector<double> {(y1[0]+y2[0])/2.0 - (x1[0]+x2[0])/2.0,(y1[1]+y2[1])/2.0 - (x1[1]+x2[1])/2.0,
                              (y1[2]+y2[2])/2.0 - (x1[2]+x2[2])/2.0};
}

//if dotproduct is > 0 then return 0 - i.e. not visible.
int face_visible (std::vector<double> x, std::vector<double> y){
  std::cout << dot(x,y)<<std::endl;
  return (dot(x,y)<0)? 1:0 ;
}

template<typename T>
void transpose(std::vector<std::vector<T>> matrix) {



}

void rotate(std::vector<double> &vector, std::vector<double> axis, double angle) {
  std::vector< std::vector<double> > R (3, std::vector<double>(3,0.));
  R[0][0] = cos(angle)+(axis[0]*axis[0])*(1.0-cos(angle));
  R[0][1] = axis[0]*axis[1]*(1.0-cos(angle) - axis[2]*sin(angle));
  R[0][2] = axis[0]*axis[2]*(1.0-cos(angle) - axis[1]*sin(angle));
  R[1][0] = axis[0]*axis[1]*(1.0-cos(angle) + axis[2]*sin(angle));
  R[1][1] = cos(angle)+(axis[1]*axis[1])*(1.0-cos(angle));
  R[1][2] = axis[1]*axis[2]*(1.0-cos(angle) - axis[0]*sin(angle));
  R[2][0] = axis[2]*axis[0]*(1.0-cos(angle) - axis[1]*sin(angle));
  R[2][1] = axis[2]*axis[1]*(1.0-cos(angle) + axis[0]*sin(angle));
  R[2][2] = cos(angle)+(axis[2]*axis[2])*(1.0-cos(angle));

  vector = R*vector;

  return;

}


int nearest_point2point(std::vector<double> point, std::vector<voxel> voxels) {

  double distance = sqrt(dot(point-voxels[0].center,point-voxels[0].center));
  int voxel_index = 0;

  for (int i=1; i<voxels.size(); i++){
    if(sqrt(dot(point-voxels[i].center,point-voxels[i].center))<distance){

      distance = sqrt(dot(point-voxels[i].center,point-voxels[i].center));
      voxel_index = i;

    }

  }
  return voxel_index;
}
