//Mesh Class implementation
#include <limits>
#include "Mesh.h"




////For the string input
//template <typename MyType>
//std::vector<MyType> String_to_Vector (std::string line){
//  std::istringstream iss(line);
//  std::vector<MyType> numbers;
//  MyType T;
//  while (iss >> T ){
//    numbers.push_back( T );
//  }
//  return numbers;
//}

// Calculate the area of a triangular face
inline double areaFace(Mesh::Face face){
  std::vector<double> a,b;
  std::vector<double> cross (3,0.);
  a = face.vertices[1]-face.vertices[0];
  b = face.vertices[2]-face.vertices[0];// 0 and 1 are an arbitrary choice of sides
  cross[0] = a[1]*b[2]-a[2]*b[1];
  cross[1] = a[0]*b[2]-a[2]*b[0];
  cross[2] = a[0]*b[1]-a[1]*b[0];

  return std::sqrt(dot(cross,cross))/2.0 ;
}

inline std::vector<double> cross_product(std::vector<double> x, std::vector<double> y){
  std::vector<double> result (x.size(),0.);
  if(x.size() == y.size()){
    for(int i=0; i<x.size(); i++){
      for(int j=0; j<y.size(); j++){
        result[i] = ((i+1)%x.size()<(i+2)%y.size()?1:-1)*x[(i+1)%x.size()]*y[(i+2)%y.size()]
                    + ((i+1)%x.size()>(i+2)%y.size()?1:-1)*y[(i+1)%x.size()]*x[(i+2)%y.size()];
      }
    }

  }else{throw("Cross_product: Arrays must be of the same size!");}
  return result;
}

// Custom print of a vector
template <typename T>
void vector_print(std::vector<T> x){
  for(int i=0; i<x.size(); i++){
      std::cout << x[i] << " ";
  }
  std::cout << std::endl;
  return;
}
//************************ Mesh Class implementation ************************//

// Mesh init implementation
Mesh::Mesh (std::string name, bool self_vf, bool opaque, double const_power_per_face){
  this->name = name;
  this->opaque = opaque;
  this->const_power_per_face = const_power_per_face;
  this->self_vf=self_vf;
}

/* Create the meshes for the simulation from a file produced via a python script acting
 * on a mesh from Blender. (Include this script for ease-of-use)
 *
 */
                //Basko values are for the scaling model.

void Mesh::CreateFromBlender(std::string filename, voxel_grid& Voxel_grid, double albedo , double basko_k , double basko_a , double basko_b){
  // Vectors to store the data for the meshes
  std::vector< std::vector<double> > vertices;
  std::vector< std::vector<int> > faceIndexes;
  std::vector< std::vector<double> > normals;
  // Variables required for reading from file
  std::string line;
  std::ifstream myfile (filename);
  bool isFirstLine = true;

  double min_area=std::numeric_limits <double> ::  max ();

  // Open the file and look through lines
  //Edit the following to change the desired mesh file format
  if(myfile.is_open()){
    while ( getline (myfile, line)){
      if(isFirstLine){
        //Convert string to vector or ints
        std::vector<int> numbers = String_to_Vector <int> (line);
        numberVertices = numbers[0];
        numberFaces = numbers[1];
        isFirstLine = false;

      }
      if(line=="v"){
        for(int i=0; i<numberVertices; i++){
          getline (myfile, line);
          vertices.push_back (std::vector <double>());
          for(int j=0; j<3; j++){
            vertices[i].push_back(String_to_Vector <double> (line)[j]);
          }

        }
      }
      if(line=="f"){
        for(int i=0; i<numberFaces; i++){
          getline (myfile, line);
          faceIndexes.push_back (std::vector <int>());
          for(int j=0; j<3; j++){
            faceIndexes[i].push_back(String_to_Vector <int> (line)[j]);
          }
        }
      }
      if(line=="n"){
        for(int i=0; i<numberFaces; i++){
          getline (myfile, line);
          normals.push_back (std::vector <double>());
          for(int j=0; j<3; j++){
            normals[i].push_back(String_to_Vector <double> (line)[j]);
          }
        }
      }

    }
    myfile.close();
  }
  else { std::cout << "Could not open file" << std::endl;}

  //Create the Face object

  for(int i=0; i<numberFaces; i++){

    Face F;
    F.nRayHits = 0;
    F.power = 0;
    F.abs_energy_flux = 0;
    F.previous_power = 0;
    F.source = this->const_power_per_face;
    F.albedo = albedo; //example value
    F.laser_eff = 1.0;

    //Basko scaling parameters
    F.a = basko_a;
    F.b = basko_b;
    F.k = basko_k;

    F.normal = normals[i];

    for(int j=0; j<faceIndexes[i].size(); j++){
      F.vertices.push_back( vertices[faceIndexes[i][j]] );
    }
    F.area = areaFace(F);
    if(F.area < min_area){min_area=F.area;} //Set the min area (used later for raytracing)
    F.areaCut = F.area;
    this->faces.push_back( F );

      //Assign as member to a voxel (- Voxels are use to speed up the ray-tracing later)
    //Center of face used to see which voxel center is closer
      std::vector<double> centroid {(this->faces[i].vertices[0][0]+this->faces[i].vertices[1][0]+this->faces[i].vertices[2][0])/3.0,
                          (this->faces[i].vertices[0][1]+this->faces[i].vertices[1][1]+this->faces[i].vertices[2][1])/3.0,
                          (this->faces[i].vertices[0][2]+this->faces[i].vertices[1][2]+this->faces[i].vertices[2][2])/3.0};
      //The voxel::member struct used to add a index, meshname pair to the voxel
      voxel::member current_face;
      current_face.face_index=i;
      current_face.mesh_name=this->name;

      //Use the nearest_point2point function in myvectorstuff to find the index of the closed voxel and push
      //the member created above onto the voxel.
      Voxel_grid.grid[nearest_point2point(centroid, Voxel_grid.grid)].content.push_back(current_face);

  }
  //Set the min area for the object
  this->min_area = min_area;
}

/*
 * The main solver for computing the VF between two surfaces
 * Employs the contour integral method over the double area integral method for increased performance
 * and accuracy.
 */
//Employs minimum distance min_dist to account for log(0)
void Mesh::viewfactor(Mesh &mesh){
  //Create the vfObject in this->
  vfObject vfObj1(this->numberFaces,mesh.numberFaces,mesh.name);
  vfObject vfObj2(mesh.numberFaces,this->numberFaces,this->name);
  //Constants
  const double pi = 3.14159265358979323846;
  const double min_dist = 1e-9;
  //Calculation variables
  std::vector< std::vector<double> > g1(
    3,
    std::vector<double> (3));//the side lengths
  std::vector<double> g2(3,0.);//the side lengths
  double m; //mid-point distance
  double sum;

  //Self-shadowing variables
  std::vector<double> x;//Vector for self_shadowing inline calculation
  std::vector< std::vector<double> > dummyverts (3,
    std::vector<double>(3,0.));
  std::vector<int> rs (3,0);//Resegmenting array
  std::vector<double> r1 (3,0.);
  std::vector<double> p1(3, 0.), p2(3, 0.);
  // Faces can block their own vision of other faces, this 'cut area' accounts for this.
  double areaCut;

    // Main face loop
    for(unsigned int i=0; i<this->numberFaces; i++){
      for(unsigned int j=0; j<mesh.numberFaces; j++){
        if(dot(this->faces[i].normal,mesh.faces[j].normal)<=0){
          //Default case
          dummyverts = mesh.faces[j].vertices;
          //AreaCut is equal to the entire face area until cut (obs)
          areaCut = mesh.faces[j].area;

          // We need to find which face of the two being considered
          // in this step of the calculation is cut by the other
          bool mesh1_cut = false;
          bool mesh2_cut = false;

          //The 'g's are the vectors pointing around the face (they form the contour)
          g1[0]= this->faces[i].vertices[1]-this->faces[i].vertices[0];
          g1[1]= this->faces[i].vertices[2]-this->faces[i].vertices[1];
          g1[2]= this->faces[i].vertices[0]-this->faces[i].vertices[2];
          //Default case ^^

          //Self-shadowing calculation for Mesh1//
          // r1 is the centroid of the face in question
          r1 = (this->faces[i].vertices[0]+
                this->faces[i].vertices[1]+
                this->faces[i].vertices[2])/3.0;
          //rs is the vector of 1s or 0s which keep track of how many faces are invisible
          //0 == visible, 1 == invisible
          rs[0]=0; rs[1]=0; rs[2]=0;
          //Are any vertices on face 2 below the face 1 horizon
          for(unsigned int l=0; l<mesh.faces[j].vertices.size(); l++){
            x = mesh.faces[j].vertices[l] - r1;
            if(dot(this->faces[i].normal,x)<0){ rs[l] = 1; mesh2_cut=true;
              //Contour segment vectors of the uncut face
              g1[0]= this->faces[i].vertices[1]-this->faces[i].vertices[0];
              g1[1]= this->faces[i].vertices[2]-this->faces[i].vertices[1];
              g1[2]= this->faces[i].vertices[0]-this->faces[i].vertices[2];

              //Dummy vertices used during the cutting of the other face
              dummyverts = mesh.faces[j].vertices;
            }
          }

          //Self-shadowing calculation for Mesh2//
          // r1 is the centroid of the face in question
          r1 = (mesh.faces[j].vertices[0]+
                mesh.faces[j].vertices[1]+
                mesh.faces[j].vertices[2])/3.0;
          //Are any vertices on face 1 below the face 2 horizon
          for(unsigned int l=0; l<this->faces[i].vertices.size(); l++){
            x = this->faces[i].vertices[l] - r1;
            if(dot(mesh.faces[j].normal,x)<0){ rs[l] = 1; mesh1_cut=true;
              //Contour segment vectors of the uncut face
              g1[0]= mesh.faces[j].vertices[1]-mesh.faces[j].vertices[0];
              g1[1]= mesh.faces[j].vertices[2]-mesh.faces[j].vertices[1];
              g1[2]= mesh.faces[j].vertices[0]-mesh.faces[j].vertices[2];

              //Dummy vertices used during the cutting of the other face
              dummyverts = this->faces[i].vertices;
            }
          }

          //Deals with two parallel planes
          if(mesh1_cut&&mesh2_cut){mesh1_cut=false;mesh2_cut=false;
              g1[0]= this->faces[i].vertices[1]-this->faces[i].vertices[0];
              g1[1]= this->faces[i].vertices[2]-this->faces[i].vertices[1];
              g1[2]= this->faces[i].vertices[0]-this->faces[i].vertices[2];

              //Dummy vertices used during the cutting of the other face
              dummyverts = mesh.faces[j].vertices;}

          //End of self-shadowing calculation

          //Apply self_shadowing by splitting segments
          //One vertex below face
          if((mesh2_cut)&&((rs[0]+rs[1]+rs[2])==1)) {
            for (int a = 0; a < 3; a++) {
              if (rs[a] == 1) {

                r1 = (this->faces[i].vertices[0]+
                      this->faces[i].vertices[1]+
                      this->faces[i].vertices[2])/3.0;

                p1 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                      dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                          this->faces[i].normal)) *
                                     (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

                p2 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                      dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a],
                                          this->faces[i].normal)) *
                                     (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a]);

                areaCut = mesh.faces[j].area - std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                             cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

                dummyverts.erase(dummyverts.begin() + a);//remove vertices
                dummyverts.insert(dummyverts.begin() + a, p2);
                dummyverts.insert(dummyverts.begin() + a, p1);

              }
            }
          }
          //Two vertices below face
          if((mesh2_cut)&&((rs[0]+rs[1]+rs[2])==2)){
            for(int a=0; a<3; a++){
              if(rs[a]==0){

                r1 = (this->faces[i].vertices[0]+
                      this->faces[i].vertices[1]+
                      this->faces[i].vertices[2])/3.0;

                p1 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                      dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                          this->faces[i].normal)) *
                                     (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

                p2 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                      dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a],
                                          this->faces[i].normal)) *
                                     (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a]);

                areaCut = std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                        cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

                dummyverts.erase(dummyverts.begin() + (a+1)%dummyverts.size());//remove vertices
                dummyverts.insert(dummyverts.begin() + (a+1)%dummyverts.size(), p1);
                dummyverts.erase(dummyverts.begin() +
                                         ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()));
                dummyverts.insert(dummyverts.begin() +
                                         ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()),p2);

              }
            }

          }

          //Apply self_shadowing by splitting segments
          //One vertex below face
          if((mesh1_cut)&&((rs[0]+rs[1]+rs[2])==1)) {
            for (int a = 0; a < 3; a++) {
              if (rs[a] == 1) {

                r1 = (mesh.faces[j].vertices[0]+
                      mesh.faces[j].vertices[1]+
                      mesh.faces[j].vertices[2])/3.0;

                p1 = dummyverts[a] + (dot(r1 - dummyverts[a], mesh.faces[j].normal) /
                                      dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                          mesh.faces[j].normal)) *
                                     (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

                p2 = dummyverts[a] + (dot(r1 - dummyverts[a], mesh.faces[j].normal) /
                                      dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a],
                                          mesh.faces[j].normal)) *
                                     (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a]);

                areaCut = this->faces[i].area - std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                                             cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

                dummyverts.erase(dummyverts.begin() + a);//remove vertices
                dummyverts.insert(dummyverts.begin() + a, p2);
                dummyverts.insert(dummyverts.begin() + a, p1);

              }
            }
          }
          //Two vertices below face
          if((mesh1_cut)&&((rs[0]+rs[1]+rs[2])==2)){
            for(int a=0; a<3; a++){
              if(rs[a]==0){

                r1 = (mesh.faces[j].vertices[0]+
                      mesh.faces[j].vertices[1]+
                      mesh.faces[j].vertices[2])/3.0;

                p1 = dummyverts[a] + (dot(r1 - dummyverts[a], mesh.faces[j].normal) /
                                      dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                          mesh.faces[j].normal)) *
                                     (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

                p2 = dummyverts[a] + (dot(r1 - dummyverts[a], mesh.faces[j].normal) /
                                      dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a],
                                          mesh.faces[j].normal)) *
                                     (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a]);

                //check that the meshes are being cut correctly
                areaCut = std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                        cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

                dummyverts.erase(dummyverts.begin() + (a+1)%dummyverts.size());//remove vertices
                dummyverts.insert(dummyverts.begin() + (a+1)%dummyverts.size(), p1);
                dummyverts.erase(dummyverts.begin() +
                                 ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()));
                dummyverts.insert(dummyverts.begin() +
                                  ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()),p2);

              }
            }

          }

          //SOLVER NUMVER 1
          //If the face is not completely self shadowed
          if(((rs[0]+rs[1]+rs[2])!=3)&&!mesh1_cut) {
            //Viewfactor solver
            for (unsigned int k = 0; k < 3; k++) {
              for (unsigned int l = 0; l < dummyverts.size(); l++) {
                g2 = dummyverts[(l + 1) % dummyverts.size()] - dummyverts[l];

                //Viewfactor calculation//
                //m: Distance between midpoints of segment vectors (min dist helps with log(0) errors)
                m = min_dist+ magnitude(midpoint(this->faces[i].vertices[(k + 1) % 3], this->faces[i].vertices[k],
                                       dummyverts[(l + 1) % dummyverts.size()], dummyverts[l]));
                sum += dot(g1[k], g2) * std::log(m);
              }
            }

            vfObj1.F[i][j] = fabs(sum) / (2.0 * pi * this->faces[i].area);
            vfObj2.F[j][i] = fabs(sum)  / (2.0 * pi * mesh.faces[j].area); //Area cut is the area of the split face that's visible

            sum = 0;
          }
          // SOLVER NUMBER 2
          else if ((mesh1_cut)&&((rs[0]+rs[1]+rs[2])!=3)){
            //Viewfactor solver
            for (unsigned int k = 0; k < 3; k++) {
              for (unsigned int l = 0; l < dummyverts.size(); l++) {
                g2 = dummyverts[(l + 1) % dummyverts.size()] - dummyverts[l];

                //Viewfactor calculation//
                //m: Distance between midpoints of segment vectors (min dist helps with log(0) errors)
                m = min_dist  + magnitude(midpoint(dummyverts[(l + 1) % dummyverts.size()], dummyverts[l],
                                       mesh.faces[j].vertices[(k + 1) % 3], mesh.faces[j].vertices[k]));
                sum += dot(g1[k], g2) * std::log(m);
              }
            }
            vfObj1.F[i][j] = fabs(sum) / (2.0 * pi * this->faces[i].area);
            vfObj2.F[j][i] = fabs(sum) / (2.0 * pi * mesh.faces[j].area); //Area cut is the area of the split face that's visible

            sum = 0;
          }
          else{vfObj1.F[i][j] = 0.;//If faces are completely below the horizon
            vfObj2.F[j][i] = 0.;}

        } else{
          vfObj1.F[i][j] = 0.;//If faces are invisible F=0
          vfObj2.F[j][i] = 0.;

        }
      }
    }
    //Push the viewfactor matrix onto the database
    this->vf_database.content.push_back( vfObj1 );
    mesh.vf_database.content.push_back( vfObj2 );

}
//Employs minimum distance min_dist to account for log(0)
void Mesh::selfviewfactor(){
  //Create the vfObject in this->
  vfObject vfObj1(this->numberFaces,this->numberFaces,this->name);
  //Constants
  const double pi = 3.14159265358979323846;
  const double min_dist = 1e-9;
  //Calculation variables
  std::vector< std::vector<double> > g1(
          3,
          std::vector<double> (3));//the side lengths
  std::vector<double> g2(3,0.);//the side lengths
  double m; //mid-point distance
  double sum;
  //Self-shadowing variables
  std::vector<double> x;//Vector for self_shadowing inline calculation

  std::vector< std::vector<double> > dummyverts (3,
                                                 std::vector<double>(3,0.));
  std::vector<int> rs (3,0);//Resegmenting array
  std::vector<double> r1 (3,0.);
  std::vector<double> p1(3, 0.), p2(3, 0.);;
  double areaCut; //For when the face is cut
  bool mesh1_cut = false;
  bool mesh2_cut = false; //Keep track of what face is cut

  for(unsigned int i=0; i<this->numberFaces; i++){
    for(unsigned int j=i; j<this->numberFaces; j++){
      //NNOTE::Try this as a multiplication factor to see if it is faster.
      if(dot(this->faces[i].normal,this->faces[j].normal)<=0){
        //Default case
        dummyverts = this->faces[j].vertices;
        areaCut = this->faces[j].area;//(i+j)%this->numberFacesUst to be sure

        bool mesh1_cut = false;
        bool mesh2_cut = false;

        g1[0]= this->faces[i].vertices[1]-this->faces[i].vertices[0];
        g1[1]= this->faces[i].vertices[2]-this->faces[i].vertices[1];
        g1[2]= this->faces[i].vertices[0]-this->faces[i].vertices[2];
        //Default case ^^

        //Self-shadowing calculation for Mesh1//
        // r1 is the centroid of the face in question
        r1 = (this->faces[i].vertices[0]+
              this->faces[i].vertices[1]+
              this->faces[i].vertices[2])/3.0;
        //rs is the vector of 1s or 0s which keep track of how many faces are invisible
        //0 == visible, 1 == invisible
        rs[0]=0; rs[1]=0; rs[2]=0;
        //Are any vertices on face 2 below the face 1 horizon
        for(unsigned int l=0; l<this->faces[j].vertices.size(); l++){
          x = this->faces[j].vertices[l] - r1;
          if(dot(this->faces[i].normal,x)<=0){ rs[l] = 1; mesh2_cut=true;
            //Contour segment vectors of the uncut face
            g1[0]= this->faces[i].vertices[1]-this->faces[i].vertices[0];
            g1[1]= this->faces[i].vertices[2]-this->faces[i].vertices[1];
            g1[2]= this->faces[i].vertices[0]-this->faces[i].vertices[2];

            //Dummy vertices used during the cutting of the other face
            dummyverts = this->faces[j].vertices;
          }
        }

        //Self-shadowing calculation for Mesh2//
        // r1 is the centroid of the face in question
        r1 = (this->faces[j].vertices[0]+
              this->faces[j].vertices[1]+
              this->faces[j].vertices[2])/3.0;
        //Are any vertices on face 1 below the face 2 horizon
        for(unsigned int l=0; l<this->faces[i].vertices.size(); l++){
          x = this->faces[i].vertices[l] - r1;
          if(dot(this->faces[j].normal,x)<=0){ rs[l] = 1; mesh1_cut=true;
            //Contour segment vectors of the uncut face
            g1[0]= this->faces[j].vertices[1]-this->faces[j].vertices[0];
            g1[1]= this->faces[j].vertices[2]-this->faces[j].vertices[1];
            g1[2]= this->faces[j].vertices[0]-this->faces[j].vertices[2];

            //Dummy vertices used during the cutting of the other face
            dummyverts = this->faces[i].vertices;
          }
        }

        //Deals with two parallel planes
        if(mesh1_cut&&mesh2_cut){mesh1_cut=false;mesh2_cut=false;
            g1[0]= this->faces[i].vertices[1]-this->faces[i].vertices[0];
            g1[1]= this->faces[i].vertices[2]-this->faces[i].vertices[1];
            g1[2]= this->faces[i].vertices[0]-this->faces[i].vertices[2];

            //Dummy vertices used during the cutting of the other face
            dummyverts = this->faces[j].vertices;}

        //End of self-shadowing calculation

        //Apply self_shadowing by splitting segments
        //One vertex below face
        if((mesh2_cut)&&((rs[0]+rs[1]+rs[2])==1)) {
          for (int a = 0; a < 3; a++) {
            if (rs[a] == 1) {

              r1 = (this->faces[i].vertices[0]+
                    this->faces[i].vertices[1]+
                    this->faces[i].vertices[2])/3.0;

              p1 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                    dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                        this->faces[i].normal)) *
                                   (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

              p2 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                    dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a],
                                        this->faces[i].normal)) *
                                   (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a]);

              areaCut = this->faces[j].area - std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                                           cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

              dummyverts.erase(dummyverts.begin() + a);//remove vertices
              dummyverts.insert(dummyverts.begin() + a, p2);
              dummyverts.insert(dummyverts.begin() + a, p1);

            }
          }
        }
        //Two vertices below face
        if((mesh2_cut)&&((rs[0]+rs[1]+rs[2])==2)){
          for(int a=0; a<3; a++){
            if(rs[a]==0){

              r1 = (this->faces[i].vertices[0]+
                    this->faces[i].vertices[1]+
                    this->faces[i].vertices[2])/3.0;

              p1 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                    dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                        this->faces[i].normal)) *
                                   (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

              p2 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[i].normal) /
                                    dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a],
                                        this->faces[i].normal)) *
                                   (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a]);

              areaCut = std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                      cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

              dummyverts.erase(dummyverts.begin() + (a+1)%dummyverts.size());//remove vertices
              dummyverts.insert(dummyverts.begin() + (a+1)%dummyverts.size(), p1);
              dummyverts.erase(dummyverts.begin() +
                               ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()));
              dummyverts.insert(dummyverts.begin() +
                                ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()),p2);

            }
          }

        }

        //Apply self_shadowing by splitting segments
        //One vertex below face
        if((mesh1_cut)&&((rs[0]+rs[1]+rs[2])==1)) {
          for (int a = 0; a < 3; a++) {
            if (rs[a] == 1) {

              r1 = (this->faces[j].vertices[0]+
                      this->faces[j].vertices[1]+
                      this->faces[j].vertices[2])/3.0;

              p1 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[j].normal) /
                                    dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                        this->faces[j].normal)) *
                                   (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

              p2 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[j].normal) /
                                    dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a],
                                        this->faces[j].normal)) *
                                   (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size())] - dummyverts[a]);

              areaCut = this->faces[i].area - std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                                            cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

              dummyverts.erase(dummyverts.begin() + a);//remove vertices
              dummyverts.insert(dummyverts.begin() + a, p2);
              dummyverts.insert(dummyverts.begin() + a, p1);

            }
          }
        }
        //Two vertices below face
        if((mesh1_cut)&&((rs[0]+rs[1]+rs[2])==2)){
          for(int a=0; a<3; a++){
            if(rs[a]==0){

              r1 = (this->faces[j].vertices[0]+
                      this->faces[j].vertices[1]+
                    this->faces[j].vertices[2])/3.0;

              p1 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[j].normal) /
                                    dot(dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a],
                                        this->faces[j].normal)) *
                                   (dummyverts[(a + 1) % dummyverts.size()] - dummyverts[a]);

              p2 = dummyverts[a] + (dot(r1 - dummyverts[a], this->faces[j].normal) /
                                    dot(dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a],
                                        this->faces[j].normal)) *
                                   (dummyverts[((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()) ] - dummyverts[a]);

              //check that the meshes are being cut correctly
              areaCut = std::sqrt(dot(cross_product(p1-dummyverts[a],p2-dummyverts[a]),
                                      cross_product(p1-dummyverts[a],p2-dummyverts[a])))/2.0;

              dummyverts.erase(dummyverts.begin() + (a+1)%dummyverts.size());//remove vertices
              dummyverts.insert(dummyverts.begin() + (a+1)%dummyverts.size(), p1);
              dummyverts.erase(dummyverts.begin() +
                               ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()));
              dummyverts.insert(dummyverts.begin() +
                                ((a - 1)<0? dummyverts.size()-1:((a-1))% dummyverts.size()),p2);

            }
          }

        }

        //SOLVER NUMVER 1
        //If the face is not completely self shadowed
        if(((rs[0]+rs[1]+rs[2])!=3)&&!mesh1_cut) {
          //Viewfactor solver
          for (unsigned int k = 0; k < 3; k++) {
            for (unsigned int l = 0; l < dummyverts.size(); l++) {
              g2 = dummyverts[(l + 1) % dummyverts.size()] - dummyverts[l];

              //Viewfactor calculation//
              //m: Distance between midpoints of segment vectors (min dist helps with log(0) errors)
              m = min_dist + magnitude(midpoint(this->faces[i].vertices[(k + 1) % 3], this->faces[i].vertices[k],
                                     dummyverts[(l + 1) % dummyverts.size()], dummyverts[l]));
              sum += dot(g1[k], g2) * std::log(m);
            }
          }

          vfObj1.F[i][j] = fabs(sum) / (2.0 * pi * this->faces[i].area);
          vfObj1.F[j][i] = fabs(sum) / (2.0 * pi * areaCut);
          sum = 0;
        }
          // SOLVER NUMBER 2
        else if ((mesh1_cut)&&((rs[0]+rs[1]+rs[2])!=3)){
          //Viewfactor solver
          for (unsigned int k = 0; k < 3; k++) {
            for (unsigned int l = 0; l < dummyverts.size(); l++) {
              g2 = dummyverts[(l + 1) % dummyverts.size()] - dummyverts[l];

              //Viewfactor calculation//
              //m: Distance between midpoints of segment vectors (min dist helps with log(0) errors)
              m = min_dist + magnitude(midpoint(dummyverts[(l + 1) % dummyverts.size()], dummyverts[l],
                                     this->faces[j].vertices[(k + 1) % 3], this->faces[j].vertices[k]));
              sum += dot(g1[k], g2) * std::log(m);
            }
          }
          vfObj1.F[i][j] = fabs(sum) / (2.0 * pi * areaCut);
          vfObj1.F[j][i] = fabs(sum) / (2.0 * pi * this->faces[j].area);

          sum = 0;
        }
        else{vfObj1.F[i][j] = 0.;vfObj1.F[j][i]=0.;//If faces are completely below the horizon
          }
          if(std::isnan(vfObj1.F[i][j])){std::cout<<"NaN in vf calc! "<<i<<" "<<j<<" "<<this->faces[i].area<<"\n"; exit(0);}
          //std::cout << "F "<< vfObj1.F[i][j]<<std::endl;

      } else{
        vfObj1.F[i][j] = 0.; vfObj1.F[j][i] = 0.;//If faces are invisible F=0

      }
    }
  }
  //Push the viewfactor matrix onto the database
  this->vf_database.content.push_back( vfObj1 );

}
//Finish doing blocking ...
void Mesh::block_with(Mesh &target_mesh, std::vector< Mesh > blocking_meshes) {
  double tolerance;
  const int nSegments = 10;
  double alpha, beta;
  std::vector<double> d (3, 0.);
  std::vector<double> v (3, 0.);
  std::vector<double> r (3, 0.);
  std::vector<double> p (3, 0.);
  double t;

    for (int i = 0; i < this->faces.size(); ++i) {
      for (int j = 0; j<target_mesh.faces.size(); ++j){
        for (int n = 0; n < blocking_meshes.size(); ++n) {
          for (int k = 0; k < blocking_meshes[n].faces.size(); ++k) {

              tolerance = sqrt(blocking_meshes[n].faces[j].area)/(nSegments);

              r = {(this->faces[i].vertices[0]+this->faces[i].vertices[1]+this->faces[i].vertices[2])/3.0};

              v = {(target_mesh.faces[k].vertices[0]+
                   target_mesh.faces[k].vertices[1]+
                    target_mesh.faces[k].vertices[2])/3.0 - r};

              if (dot(v, blocking_meshes[n].faces[k].normal) < 0) {
                t = dot((blocking_meshes[n].faces[k].vertices[0] +
                        blocking_meshes[n].faces[k].vertices[1] +
                        blocking_meshes[n].faces[k].vertices[2]) / 3.0 - r,
                        blocking_meshes[n].faces[k].normal) / dot(v, blocking_meshes[n].faces[k].normal);

                p = t * v + r;

              }

            for(int a = 0; a<=nSegments; a++) {
              alpha = 1.0 * a / (1.0 * nSegments);
              for(int b = 0; b<=nSegments-a; b++) {
                beta = 1.0 * b / (1.0 * nSegments);
                d = (alpha * ( blocking_meshes[n].faces[k].vertices[1] - blocking_meshes[n].faces[k].vertices[0]) +
                     beta * ( blocking_meshes[n].faces[k].vertices[2] - blocking_meshes[n].faces[k].vertices[0]))
                    -(p- blocking_meshes[n].faces[k].vertices[0]);

                if (magnitude(d) < tolerance) {

                  this->vf_database[target_mesh.name].F[i][j]=0.;
                  target_mesh.vf_database[this->name].F[j][i]=0.;

                  goto ray_hit;
                }
              }
            }

            for(int b = 0; b<=nSegments; b++) {
              beta = 1.0 * b / (1.0 * nSegments);
              for(int a = 0; a<=nSegments-a; a++) {
                alpha = 1.0 * b / (1.0 * nSegments);
                d = (alpha * (blocking_meshes[n].faces[k].vertices[1] - blocking_meshes[n].faces[k].vertices[0]) +
                     beta * (blocking_meshes[n].faces[k].vertices[2] - blocking_meshes[n].faces[k].vertices[0]))
                    -(p-blocking_meshes[n].faces[k].vertices[0]);

                if (magnitude(d) < tolerance) {

                  this->vf_database[target_mesh.name].F[i][j]=0.;
                  target_mesh.vf_database[this->name].F[j][i]=0.;

                  goto ray_hit;
                }
              }
            }

          }
      }
      ray_hit:
        continue;
    }
  }
}

void Mesh::print_vf(std::string meshname){
  for(int i=0; i<this->vf_database.content.size(); i++){
    if(this->vf_database.content[i].otherMesh == meshname){
      for(int j=0; j<this->vf_database.content[i].F.size(); j++){
        for(int k=0; k<this->vf_database.content[i].F[j].size(); k++){
          std::cout << "VF: " << this->vf_database.content[i].F[j][k] << std::endl;
        }
      }
    }
  }
}

void Mesh::print_totalvf(Mesh mesh){
  double lsum = 0.;
  double ksum = 0.;
  double area1 = 0.;
  double area2 = 0.;
  //Firstly, get the area of the whole of the
  //second mesh (A_j)
  for(int i=0; i<mesh.faces.size(); i++){
    area2 += mesh.faces[i].area;
  }
  for(int i=0; i<this->faces.size(); i++) {
    area1 += this->faces[i].area;
  }

  //Area of the first mesh is calculated in the loop below
  for(int n=0; n<this->vf_database.content.size(); n++){
    if(this->vf_database.content[n].otherMesh == mesh.name){
      for(int k=0; k<this->vf_database.content[n].F.size(); k++){
        for(int l=0; l<this->vf_database.content[n].F[k].size(); l++){
            if(this->vf_database.content[n].F[k][l]<0){ std::cout << "VF<0"<<std::endl; }
          lsum += this->vf_database.content[n].F[k][l];
            if(std::isnan(lsum)){std::cout<<"NaN in vf database\n"; exit(0);}
        }
        if(this->faces[k].area<0){ std::cout << "area<0" <<std::endl; }
        ksum += lsum*(this->faces[k].area);
        //std::cout<< ksum<<std::endl;
        lsum = 0.;

      }
    } else{
      //std::cout << "There is no VF database entry."<<std::endl;
    }
  }
  std::cout << this->name <<" from " << mesh.name
            <<" "<<ksum/(area1) << std::fixed << std::setprecision(8) << std::endl;
}

double Mesh::total_area() {
  double total_area = 0;
  for(int i=0; i<this->faces.size(); i++){
    total_area += this->faces[i].area;
  }
  return total_area;
}

/* Finish later - for use with Voxel method (if employed)
std::vector<int> Mesh::findinVolume(std::vector<double> origin[3], double width) {
  std::vector<int> result;
  for(int i=0; i<this->faces.size(); ++i){
    for(int j=0; j<this->faces[i].vertices.size(); ++j){
      if(this->faces[i].vertices[j][0]<=origin[0]+width/2.0)
    }
  }
  return std::vector<int>();
}
*/

Mesh::vfObject::vfObject(int n, int m, std::string otherMesh, double default_values) {
  this->otherMesh = otherMesh;
  F = std::vector< std::vector<double> > (n,
          std::vector<double> (m,0.));

}



void saveMeshData(Mesh mesh, int timestep, double realtime ) {
  std::ofstream save ("./pwsrc_"+mesh.name+"_"+std::to_string(timestep)+".vtk");
  save << "# vtk DataFile Version 2.0"<< std::endl;
  save << mesh.name +"\nASCII\nDATASET POLYDATA"<<std::endl;
  save << "POINTS "<< 3*mesh.numberFaces << " float"<<std::endl;
  for (int i =0; i<mesh.faces.size(); i++) {
    for (int j=0; j<3; j++) {
      save << mesh.faces[i].vertices[j][0] << " " <<
           mesh.faces[i].vertices[j][1]
           << " " << mesh.faces[i].vertices[j][2] << std::endl;
    }
  }
  save << "POLYGONS " << mesh.numberFaces << " " << (1+3)*mesh.numberFaces <<std::endl;
  for (int i=0; i<mesh.numberFaces; i++){
    save << "3 "<< 3*i << " " <<3*i + 1 << " " << 3*i + 2 << std::endl;
  }
  save << "CELL_DATA " << mesh.numberFaces << std::endl;
  save << "SCALARS power float 1\nLOOKUP_TABLE default\n";
  for (int i=0; i<mesh.numberFaces; i++){
    save << mesh.faces[i].power*(mesh.faces[i].area) << std::endl;
  }
    save << "SCALARS intensity float 1\nLOOKUP_TABLE default\n";
    for (int i=0; i<mesh.numberFaces; i++){
        save << mesh.faces[i].power<< std::endl;
    }
    save << "SCALARS source float 1\nLOOKUP_TABLE default\n";
    for (int i=0; i<mesh.numberFaces; i++){
        save << mesh.faces[i].source<< std::endl;
    }
  save << "SCALARS abs_energy_flux float 1\nLOOKUP_TABLE default\n";
  for (int i=0; i<mesh.numberFaces; i++){
    save << mesh.faces[i].abs_energy_flux<< std::endl;
  }
  save << "SCALARS albedo float 1\nLOOKUP_TABLE default\n";
  for (int i=0; i<mesh.numberFaces; i++){
    save << mesh.faces[i].albedo << std::endl;
  }
    save << "SCALARS tr float 1\nLOOKUP_TABLE default\n";
    for (int i=0; i<mesh.numberFaces; i++){
        save << sqrt(sqrt(mesh.faces[i].power/(5.67e-8)))/(11604.5250061657) << std::endl;
    }


  save.close();
    //*******************  Save Average albedo over time  *********
    std::ofstream save2 ("./albedo_vs_t_"+mesh.name+".txt");
    double total=0.;
    save2 << mesh.name+"\n# time(s)  average albedo\n";
    for (int i = 0; i < mesh.faces.size(); i++){
      total += mesh.faces[i].albedo;
    }
    save2 << realtime <<" "<< total/mesh.faces.size() << std::endl;

    save2.close();

  //*******************  Save Average Tr over time  *********
  std::ofstream save3 ("./Tr_vs_t_"+mesh.name+".txt");
  total=0.;
  save3 << mesh.name+"\n# time(s)  average albedo\n";
  for (int i = 0; i < mesh.faces.size(); i++){
    total += sqrt(sqrt(mesh.faces[i].power/(5.67e-8)))/(11604.5250061657);
  }
  save3 << realtime <<" "<< total/mesh.faces.size() << std::endl;

  save3.close();
}

//************************ Scene Class implementation ************************//
