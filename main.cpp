/* Main code for the C++ version of the FORTRAN HOHVF
- Joe Owen
*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <time.h>
#include <iomanip> //for setprecision

#include "headers/Mesh.h"
#include "headers/raytracing.h"
#include "headers/radiosity.h"
#include "headers/voxels.h"

double ray_power_global=1;

std::vector<double> operator+(std::vector<double> x,std::vector<double> y){
  return std::vector<double> {x[0]+y[0],x[1]+y[1], x[2]+y[2]};
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

class Scene{
public:
    std::vector<Mesh> meshes;
    std::vector<Mesh> blocking_meshes;
    std::vector<Ray> rays;
    Mesh& operator[] (std::string name){
      for(int i=0; i<meshes.size(); i++){
        if(meshes[i].name==name){return meshes[i];}
      }
      for(int i=0; i<blocking_meshes.size(); i++){
          if(blocking_meshes[i].name==name){return blocking_meshes[i];}
      }
      throw "Name not in array!";
    }
};


void saveSourcedata(Mesh mesh){
        std::ofstream save("./source_" + mesh.name + ".save");

        for (int i = 0; i < mesh.faces.size(); i++) {
                save << mesh.faces[i].source <<std::endl;
            }

        save.close();

    return;
}

void loadSourcedata(Mesh mesh){
    std::ifstream load ("./source_"+ mesh.name + ".save");
    std::string line;

    for (int i=0; i<mesh.faces.size(); i++){

            getline(load, line);
            std::istringstream iss(line);
            double source;
            iss >> source;
            mesh.faces[i].source = source;

    }

    load.close();

    return;
}


void saveVFdata(Mesh mesh){
    for (int k=0; k<mesh.vf_database.content.size(); k++) {
        std::ofstream save("./vf_" + mesh.name + "_" + mesh.vf_database.content[k].otherMesh + ".save");
        std::vector<std::vector<double> > F(mesh.numberFaces, std::vector<double>(mesh.vf_database.content[k].F[0].size(), 0.));
        F = mesh.vf_database.content[k].F;
        for (int i = 0; i < mesh.vf_database.content[k].F.size(); i++) {
            for (int j=0; j < mesh.vf_database.content[k].F[i].size(); j++){
                save << F[i][j] <<std::endl;
            }


        }
        save.close();
    }
    return;
}

void loadVFdata(Mesh &mesh, Mesh othermesh){
    std::ifstream load ("./vf_"+mesh.name+"_"+othermesh.name+".save");
    Mesh::vfObject V (mesh.numberFaces, othermesh.numberFaces, othermesh.name);
    std::string line;

    for (int i=0; i<mesh.numberFaces; i++){
        for (int j=0; j<othermesh.numberFaces; j++){
            getline(load, line);
            std::istringstream iss(line);
            iss >> V.F[i][j];
        }
    }

    mesh.vf_database.content.push_back(V);
    load.close();

    return;
}

std::vector< std::vector<double> > loadTVdata(std::string filename){
    std::ifstream load ("./"+filename);
    std::vector< std::vector<double> > tv_record;

    std::string line;

    while (getline(load, line)){
        std::istringstream iss(line);
        tv_record.push_back(String_to_Vector<double>(line));

    }
    return tv_record;
}

const std::string cmd_options[] = {"-h","-m","-r","-b","-mb","-ms",""};

int main(int argc, char* argv[]) {
    //Initialise the voxel grid into main for creation within the following user input stage
    voxel_grid Voxel_grid = create_voxel_grid(std::vector<double>{0.,0.,0.},1,1,1,1);
    bool use_voxels = false;//Initialize to false until set in the user input
    std::vector< std::vector<double> > tv_record;
    bool use_tv_record = false;

  if(argc==1){std::cout<<"Help file should go here... tbc... :)\n";exit(0);}
  for (int i=1; i<argc; ++i){
    if(std::string(argv[i])=="-h"){
      std::cout<<"Help:\n Info on program goes here...\n"
                 "\n-h brings up this help text.\n"
                 "\n-m indicates that a mesh record is to follow with format [name] [path/to/file]\n"
                      "\ti.e. -m name1 /path/to/mesh1 .\n"
                 "\n-mb indicates that a record for a blocking (opaque) mesh is to follow. These meshes ARE included in view factor/radiosity calculations\n"
                 "\n-ms indicates that a record for a convex (or self-viewing) mesh is to follow.\n"
                 "\n Additional options for -m, -mb and -ms records are as follows:\n"
                      "\t--raytracing meaning rays will be able to illuminate this object.\n"
                      "\t--basko k a b - these signify that the Basko scaling is to be used for the walls,k, a and b are the parameters"
                 "\n-r indicates that a ray record is to follow with format [path/to/file]\n"
                 "\n-b indicates that a record for a blocking (opaque) mesh is to follow. These 'b-meshes' are NOT included in view factor/radiosity calculations\n"
                 "\n-s indicates that a record for a diffuse source is to follow with format [source value] [name] [path/to/file].\n"
                 "\n--albedo specifies the albedo for the name mesh (MUST BE AFTER THE NAMED MESH RECORD), --albedo [name] [albedo value]"
                 "\n--laseff specifies laser X-ray conversion efficiencies (MUST BE AFTER THE NAMED MESH RECORD) , --laseff [name] [laseff value](MUST BE AFTER THE NAMED MESH RECORD)"
                 "\n--voxels use the voxel raytracer, (MUST BE BEFORE ANY MESH RECORDS) --voxels [x_width] [y_width] [z_width] [voxels per unit length]"
                 " These s-meshes ARE included in view factor/radiosity calculations.\n"
                 "\n -tv is a time value record for the lasers"
                 "\n Info on file formats to follow...\n" ;exit(0);
    }
    if(argv[i]=="-m"){
      std::cout<<"-m indicates that mesh records follow\n";exit(0);
    }
  }



  clock_t start = clock();
  //Display #1
  std::cout << "View Factor Radiosity Solver"<< std::endl;
  std::cout << "----------------------------"<<std::endl;

  //Create the scene - implement as class in the future
  Scene scene;

  //Perhaps put all of this in a seperate header file for reading in from the command line
    //#################################################################################
  //Read in meshes
  for(int i=1; i<argc; ++i) {
        //VOXEL RECORD MUST BE READ IN FIRST - allocation to the voxel array occurs within the CreateFromBlender function
        //So it needs to exist first otherwise the "empty" one will be used.
      if(std::string(argv[i])=="--voxels"){

          //Input error checking...
          //Check for double option inputs/ empty fields
          try {
              for (auto x : cmd_options) {
                  if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
                      std::cout << "Invalid input (see help, -h)\n";
                      exit(0);
                  }
              }
          }
          catch (std::logic_error) {std::cout << "Invalid input (see help, -h)\n";exit(0);}
          //End of input error checking.

          std::istringstream A (argv[i+1]);
          double x_width;
          A>>x_width;

          std::istringstream B (argv[i+2]);
          double y_width;
          B>>y_width;

          std::istringstream C (argv[i+3]);
          double z_width;
          C>>z_width;

          std::istringstream D (argv[i+4]);
          double n_voxels;
          D>>n_voxels;

          //Before anything, create the voxel grid to speed up the raytracing and occultation steps
          //Create voxel grid
          Voxel_grid = create_voxel_grid(std::vector<double>{0.,0.,0.},x_width,y_width,z_width,n_voxels);
          use_voxels = true;

      }


    if(std::string(argv[i])=="-m") {
      //Input error checking...
      //Check for double option inputs/ empty fields
        try {
          for (auto x : cmd_options) {
            if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
              std::cout << "Invalid input (see help, -h)\n";
              exit(0);
            }
          }
        }
        catch (std::logic_error) {std::cout << "Invalid input (see help, -h)\n";exit(0);}
      //End of input error checking.

      //Options
        if(std::string(argv[i+1])=="--raytracing" and std::string(argv[i+2])=="--basko")
        {
            std::cout << "Reading in " << argv[i+6] << " from "<<argv[i+7]<<" ..." << std::endl;
            scene.meshes.push_back(Mesh(std::string(argv[i+6]), false, false, 0.));
            scene[argv[i+6]].raytracing = true;
            scene[argv[i+6]].basko_scaling = true;
            //Create mesh objects from the blender input files
            //Check for scaling option

            //Read in values

            std::istringstream A (argv[i+3]);
            double k;
            A>>k;

            std::istringstream B (argv[i+4]);
            double a;
            B>>a;

            std::istringstream C (argv[i+5]);
            double b;
            C>>b;

            scene[argv[i+6]].CreateFromBlender(argv[i+7], Voxel_grid, 1., k, a, b );

        }
        else if(std::string(argv[i+1])=="--basko")
        {
            std::cout << "Reading in " << argv[i+5] << " from "<<argv[i+6]<<" ..." << std::endl;
            scene.meshes.push_back(Mesh(std::string(argv[i+5]), false, false, 0.));
            scene[argv[i+5]].basko_scaling = true;
            //Create mesh objects from the blender input files
            //Check for scaling option

            //Read in values

            std::istringstream A (argv[i+2]);
            double k;
            A>>k;

            std::istringstream B (argv[i+3]);
            double a;
            B>>a;

            std::istringstream C (argv[i+4]);
            double b;
            C>>b;

            scene[argv[i+5]].CreateFromBlender(argv[i+6], Voxel_grid, 1., k, a, b );

        }
      else if(std::string(argv[i+1])=="--raytracing" ){

          std::cout << "Reading in " << argv[i+2] << " from "<<argv[i+3]<<" ..." << std::endl;
          scene.meshes.push_back(Mesh(std::string(argv[i+2]), false, false, 0.));
          scene[argv[i+2]].raytracing = true;
          //Create mesh objects from the blender input files
          //Check for scaling option

          scene[argv[i+2]].CreateFromBlender(argv[i+3], Voxel_grid);
        }
        else{

          std::cout << "Reading in " << argv[i+1] << " from "<<argv[i+2]<<" ..." << std::endl;
          scene.meshes.push_back(Mesh(std::string(argv[i+1]), false, false, 0.));
          //Create mesh objects from the blender input files
          scene[argv[i+1]].CreateFromBlender(argv[i+2], Voxel_grid);
        }
        //Additional options
    }
        if(std::string(argv[i])=="-ms") {
          //Input error checking...
          //Check for double option inputs/ empty fields
          try {
              for (auto x : cmd_options) {
                  if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
                      std::cout << "Invalid input (see help, -h)\n";
                      exit(0);
                  }
              }
          }
          catch (std::logic_error) {std::cout << "Invalid input (see help, -h)\n";exit(0);}
          //End of input error checking.

          //Options
            if(std::string(argv[i+1])=="--raytracing" and std::string(argv[i+2])=="--basko")
            {
                std::cout << "Reading in " << argv[i+6] << " from "<<argv[i+7]<<" ..." << std::endl;
                scene.meshes.push_back(Mesh(std::string(argv[i+6]), true, false, 0.));
                scene[argv[i+6]].raytracing = true;
                scene[argv[i+6]].basko_scaling = true;
                //Create mesh objects from the blender input files
                //Check for scaling option

                //Read in values

                std::istringstream A (argv[i+3]);
                double k;
                A>>k;

                std::istringstream B (argv[i+4]);
                double a;
                B>>a;

                std::istringstream C (argv[i+5]);
                double b;
                C>>b;

                scene[argv[i+6]].CreateFromBlender(argv[i+7], Voxel_grid, 1., k, a, b );

            }
            else if(std::string(argv[i+1])=="--basko")
            {
                std::cout << "Reading in " << argv[i+5] << " from "<<argv[i+6]<<" ..." << std::endl;
                scene.meshes.push_back(Mesh(std::string(argv[i+5]),true, false, 0.));
                scene[argv[i+5]].basko_scaling = true;
                //Create mesh objects from the blender input files
                //Check for scaling option

                //Read in values

                std::istringstream A (argv[i+2]);
                double k;
                A>>k;

                std::istringstream B (argv[i+3]);
                double a;
                B>>a;

                std::istringstream C (argv[i+4]);
                double b;
                C>>b;

                scene[argv[i+5]].CreateFromBlender(argv[i+6], Voxel_grid, 1., k, a, b );

            }
          else if(std::string(argv[i+1])=="--raytracing"){

              std::cout << "Reading in " << argv[i+2] << " from "<<argv[i+3]<<" ..." << std::endl;
              scene.meshes.push_back(Mesh(std::string(argv[i+2]), true, false, 0.));
              scene[argv[i+2]].raytracing = true;
              //Create mesh objects from the blender input files
              scene[argv[i+2]].CreateFromBlender(argv[i+3], Voxel_grid);
          }
          else{

              std::cout << "Reading in " << argv[i+1] << " from "<<argv[i+2]<<" ..." << std::endl;
              scene.meshes.push_back(Mesh(std::string(argv[i+1]), true, false, 0.));
              //Create mesh objects from the blender input files
              scene[argv[i+1]].CreateFromBlender(argv[i+2], Voxel_grid);
          }
          //Additional options
      }
      if(std::string(argv[i])=="-s") {

          //Input error checking...
          //Check for double option inputs/ empty fields
          try {
              for (auto x : cmd_options) {
                  if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
                      std::cout << "Invalid input (see help, -h)\n";
                      exit(0);
                  }
              }
          }
          catch (std::logic_error) {
              std::cout << "Invalid input (see help, -h)\n";
              exit(0);
          }
          //End of input error checking.

          std::istringstream A(argv[i + 1]);
          double source;
          A >> source;

          std::cout << "Reading in " << argv[i + 2] << " from " << argv[i + 3] << " ..." << std::endl;
          scene.meshes.push_back(Mesh(std::string(argv[i + 2]), false, false, source));
          //Create mesh objects from the blender input files
          scene[argv[i + 2]].CreateFromBlender(argv[i + 3], Voxel_grid);
      }
      if(std::string(argv[i])=="-b") {

          //Input error checking...
          //Check for double option inputs/ empty fields
          try {
              for (auto x : cmd_options) {
                  if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
                      std::cout << "Invalid input (see help, -h)\n";
                      exit(0);
                  }
              }
          }
          catch (std::logic_error) {
              std::cout << "Invalid input (see help, -h)\n";
              exit(0);
          }
          //End of input error checking.

          std::cout << "Reading in blocking mesh " << argv[i + 1] << " from " << argv[i + 2] << " ..." << std::endl;
          scene.blocking_meshes.push_back(Mesh(std::string(argv[i + 1]), false, false, 0.));
          //Create mesh objects from the blender input files
          scene[argv[i + 1]].CreateFromBlender(argv[i + 2], Voxel_grid);
      }
          //Additional options

          if(std::string(argv[i])=="--albedo"){

              //Input error checking...
              //Check for double option inputs/ empty fields
              try {
                  for (auto x : cmd_options) {
                      if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
                          std::cout << "Invalid input (see help, -h)\n";
                          exit(0);
                      }
                  }
              }
              catch (std::logic_error) {std::cout << "Invalid input (see help, -h)\n";exit(0);}
              //End of input error checking.

              std::string name = std::string(argv[i+1]);

              std::istringstream A (argv[i+2]);
              double albedo;
              A>>albedo;

              for (int i=0; i<scene[name].faces.size();i++) {
                  scene[name].faces[i].albedo = albedo;
              }
          }

      if(std::string(argv[i])=="--laseff"){

          //Input error checking...
          //Check for double option inputs/ empty fields
          try {
              for (auto x : cmd_options) {
                  if (std::string(argv[i + 1]) == x || std::string(argv[i + 2]) == x) {
                      std::cout << "Invalid input (see help, -h)\n";
                      exit(0);
                  }
              }
          }
          catch (std::logic_error) {std::cout << "Invalid input (see help, -h)\n";exit(0);}
          //End of input error checking.

          std::string name = std::string(argv[i+1]);

          std::istringstream A (argv[i+2]);
          double laseff;
          A>>laseff;

          for (int i=0; i<scene[name].faces.size();i++) {
              scene[name].faces[i].laser_eff = laseff;
          }
      }

          if(std::string(argv[i])=="-tv"){

              tv_record = loadTVdata(std::string(argv[i+1]));
              use_tv_record = true;
      }


  }
  //############################################################################################

  //Read in rays
  for(int i=1; i<argc; ++i) {
    if (std::string(argv[i]) == "-r") {

      //Input error checking...
      //Check for double option inputs/ empty fields
      try {
        for (auto x : cmd_options) {
          if (std::string(argv[i + 1]) == x) {
            std::cout << "Invalid input (see help, -h)\n";
            exit(0);
          }
        }
      }
      catch (std::logic_error) {std::cout << "Invalid input (see help, -h)\n";exit(0);}
      //End of input error checking.

      //Create Ray vector from the text input file
      scene.rays = rays_from_txt(argv[i+1], ray_power_global);

      //Normalise tv record to number of rays
      if(use_tv_record){
          for (int i=0; i<tv_record.size(); i++){
              tv_record[i][1] /= scene.rays.size();
          }
      }
    }
  }

  //Timing stop
  clock_t stop = clock();
  double elapsed = (double) (stop - start)/CLOCKS_PER_SEC;
  std::cout << "Time taken: " << elapsed << " s"<<std::endl;
  //Restart timer
  start = clock();

  //Display #2 - Mesh information
  for(int i=0; i<scene.meshes.size(); ++i) {
    std::cout << scene.meshes[i].name << ": Number of faces = " << scene.meshes[i].faces.size() << std::endl;
  }

  //Display the maximum, minimum and average faces per voxel.
  int voxMax =0;
  int voxMin = Voxel_grid.grid[0].content.size();
  double voxMean =0;

  for (int k=0; k< Voxel_grid.grid.size(); k++){

        //Initialize voxMin
      if(Voxel_grid.grid[k].content.size() < voxMin){
          voxMin = Voxel_grid.grid[k].content.size();
      }
      if(Voxel_grid.grid[k].content.size()>voxMax){
            voxMax = Voxel_grid.grid[k].content.size();
        }
      voxMean += Voxel_grid.grid[k].content.size();
  }
  voxMean /= Voxel_grid.grid.size();

  if(voxMean*Voxel_grid.grid.size() < scene.meshes[0].faces.size() + scene.meshes[1].faces.size()){
      std::cout << "Not enough voxels!!\n"; exit(0);
  }
  std::cout << "Average # faces per voxel = " << voxMean <<std::endl;
  std::cout << "Max # faces per voxel = " << voxMax <<std::endl;
  std::cout << "Min # faces per voxel = " << voxMin <<std::endl;
//########################################################################

  char cbuffer;
  std::cout<< "Would you like to load the view factors for these meshes? (y/n): ";
  std::cin >> cbuffer;
  //Load view factors option
  if(cbuffer=='y'){
      std::cout<<"Loading view factor data...\n";
      for (int i = 0; i < scene.meshes.size(); ++i) {
          for (int j = 0; j < scene.meshes.size(); ++j) {
              //Do the self vf if i==j
              if (i == j && scene.meshes[i].self_vf) { loadVFdata(scene.meshes[i], scene.meshes[i]); }
                  //Calc view factors with care taken to take advantage of reciprocity relation.
              if (i != j) { loadVFdata(scene.meshes[i], scene.meshes[j]); }
          }
      }
      std::cout<<"...Done!\n";
  }
  else {
      //View factor calculations + save
      std::cout<< "Calculating view factors...\n";
      for (int i = 0; i < scene.meshes.size(); ++i) {
          for (int j = i; j < scene.meshes.size(); ++j) {
              //Do the self vf if i==j
              if (i == j && scene.meshes[i].self_vf) { scene.meshes[i].selfviewfactor(); }
                  //Calc view factors with care taken to take advantage of reciprocity relation.
              else if (i != j) { scene.meshes[i].viewfactor(scene.meshes[j]); }
          }
      }
  }
    //End of loading / calculating viewfactors

  //Blocking from bmeshes
    if(scene.blocking_meshes.size()>0) {
        for (int i = 0; i < scene.meshes.size(); i++) {
            for (int j = 0; j < scene.meshes.size(); ++j) {
                //Do the self vf blocking if i==j
                if (i == j && scene.meshes[i].self_vf) {
                    scene.meshes[i].block_with(scene.meshes[i], scene.blocking_meshes);
                }
                    //Calc blocking with care taken to take advantage of reciprocity relation.
                else if (i != j) { scene.meshes[i].block_with(scene.meshes[j], scene.blocking_meshes); }
            }
        }
    }


    std::cout << "Displaying total view factors\n";
    //Display #3 - View Factors
    for (int i = 0; i < scene.meshes.size(); ++i) {
        for (auto vf : scene.meshes[i].vf_database.content) {
            scene.meshes[i].print_totalvf(scene[vf.otherMesh]);
        }
    }

    if (cbuffer != 'y') {
        std::cout << "Saving view factors ...\n";
        //Save view factors
        for (int i = 0; i < scene.meshes.size(); i++) {
            saveVFdata(scene.meshes[i]);
        }
    }



  //Perform the raytracing (now with voxel option)

    std::cout<< "Would you like to load the source data for the meshes? (y/n): ";
    std::cin >> cbuffer;

    if(cbuffer == 'y') {

        for(int i =0; i<scene.meshes.size(); ++i){
            if(scene.meshes[i].raytracing){
                loadSourcedata(scene.meshes[i]);
            }
        }

    }
    else {
        std::cout<<"Performing raytracing... \n";
        for (int i = 0; i < scene.meshes.size(); ++i) {
            if (scene.meshes[i].raytracing && !use_voxels) { raytrace(scene.rays, scene.meshes[i]); }
            if (scene.meshes[i].raytracing && use_voxels) {
                voxel_raytrace_nearest(scene.rays, scene.meshes[i], Voxel_grid);
            }

        }
        std::cout << "Would you like to save the source data for these meshes? (y/n): ";
        std::cin >> cbuffer;
        if(cbuffer == 'y'){
            for(int k =0; k<scene.meshes.size(); k++){
                if(scene.meshes[k].raytracing){
                    saveSourcedata(scene.meshes[k]);
                }
            }
        }
    }

  std::cout<<"radiosity\n";


    if(scene.meshes[0].basko_scaling and scene.meshes[1].basko_scaling){
        time_dependent_solveRadiosity(scene.meshes[0],scene.meshes[1], tv_record);
    }
    else{
        solveRadiosity(scene.meshes[0],scene.meshes[1]);

        for(int i =0; i<scene.meshes.size(); ++i) {
            saveMeshData(scene.meshes[i]);

        }
    }



    for(int i =0; i<scene.meshes.size(); ++i) {
        std::cout << "Area " << scene.meshes[i].name << " = " << scene.meshes[i].total_area() << std::endl;
    }

  //Timing stop
    stop = clock();

  elapsed = (double) (stop - start)/CLOCKS_PER_SEC;
  std::cout << "Time taken: " << elapsed << std::endl;

  return 0;
}
