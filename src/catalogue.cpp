// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "parameters.h"
#include "catalogue.h"



// ==========================================================
// normalize vectors in input catalogue (allow for interpretation of inner products in terms of angles)
// ==========================================================
void catalogue::normalize_vectors(){
    
    for(int i=0; i<input.obj.size();i++){
        
        double vec_a_abs = sqrt( std::inner_product(std::begin(input.obj[i].vec_a), std::end(input.obj[i].vec_a), std::begin(input.obj[i].vec_a), 0.0) );
        double vec_b_abs = sqrt( std::inner_product(std::begin(input.obj[i].vec_b), std::end(input.obj[i].vec_b), std::begin(input.obj[i].vec_b), 0.0) );
                
        for(int j=0; j<input.obj[i].vec_a.size(); j++){
            input.obj[i].vec_a[j] /= vec_a_abs;
        }
        
        for(int j=0; j<input.obj[i].vec_b.size(); j++){
            input.obj[i].vec_b[j] /= vec_b_abs;
        }
    }
}



// ==========================================================
// read catalogue from ascii file
// ==========================================================
void catalogue::read(const parameters p, const std::string filename){

    
    std::ifstream indata;
    indata.open(filename.c_str());
    

    if(!indata.is_open()){
        std::cout<<"##### ERROR [read_cat]: can't open: " << filename <<" #####"<<std::endl;
        exit (EXIT_FAILURE);
    }else{
        
        //------- loop over lines in input file -------
        std::string line;
        while (std::getline(indata, line)){
        
            std::stringstream  ss(line);
            std::string var;
            int columnNr = 0;    
            object obj_line;
            
            while(std::getline(ss,var, p.delim_in[0])){

                //check input
                if(std::isnan(stod(var))){
                    std::cerr<< "#### ERROR [cat.read]: NaN in line " << input.obj.size() << ", column " << columnNr <<" ####" <<std::endl;
                    exit(EXIT_FAILURE);
                }
                if(std::isinf(stod(var))){
                    std::cerr<< "#### ERROR [cat.read]: inf in line " << input.obj.size() << ", column " << columnNr << " ####" << std::endl;
                    exit(EXIT_FAILURE);
                }

                //check if colunm is in input column list
                auto found_pos = std::find(std::begin(p.cols_pos), std::end(p.cols_pos), columnNr);
                auto found_vec_a = std::find(std::begin(p.cols_vec_a), std::end(p.cols_vec_a), columnNr);
                auto found_vec_b = std::find(std::begin(p.cols_vec_b), std::end(p.cols_vec_b), columnNr);
                
                
                //add variable to data vectors
                if (found_pos != std::end(p.cols_pos)) { obj_line.pos.push_back(std::stod(var)); }
                if (found_vec_a != std::end(p.cols_vec_a)) { obj_line.vec_a.push_back(std::stod(var)); }
                if (found_vec_b != std::end(p.cols_vec_b)) { obj_line.vec_b.push_back(std::stod(var)); }
        
                columnNr++;
            }
            input.obj.push_back(obj_line);
        }
    }
    
}



// ==========================================================
// find minima and maxima of positions
// ==========================================================
std::vector < std::vector < double > > catalogue::limits(sample smp){

    std::vector < std::vector < double > > lim;

    //initilize min/max position with first position vector
    std::vector < double > init = smp.obj[0].pos;
    
    for(int i=0; i< init.size(); i++){
        lim.push_back({init[i],init[i]});
    }

    //search for min/max
    for(int i=0; i< smp.obj.size(); i++){
        for(int j=0; j<smp.obj[i].pos.size(); j++){
            
            if(smp.obj[i].pos[j] < lim[j][0]){
                lim[j][0] = smp.obj[i].pos[j];
            }

            if(smp.obj[i].pos[j] > lim[j][1]){
                lim[j][1] = smp.obj[i].pos[j];
            }
        }  
    }
    
    return lim;
    
    //std::cout<<pos_min[0]<<"\t"<<pos_max[0]<<std::endl;
}



//===========================================================
// convert cellID to grid coordinates
//===========================================================
//WARNING: max cells per axis = (4294967295)^(1÷3) = 1625
// void ID_to_pos(const unsigned long ID, const int DIM, int pos[3]){
//     
//     pos[2] = int(ID / DIM / DIM);
//     pos[1] = int(ID / DIM) - pos[2]*DIM;
//     pos[0] = ID - pos[1]*DIM - pos[2]*DIM*DIM;
// }



//===========================================================
// convert cellID to grid coordinates
//===========================================================
//WARNING: max cells per axis = (4294967295)^(1÷3) = 1625
// unsigned long pos_to_ID(const int DIM, const int I, const int J, const int K){    
//     return I + J*DIM + K*DIM*DIM;
// }



// ==========================================================
// obtain cell ID for object from its position
// ==========================================================
int catalogue::pos_to_ID(
    const std::vector < int > & numb_jk,
    const std::vector < double > & pos,
    const std::vector < std::vector < double > > & pos_limits,
    const std::vector < double > & Lcell){

    int jk_coords[pos.size()]={0};

    for(int i=0; i< pos.size(); i++){
        
        jk_coords[i] = int( (pos[i] - pos_limits[i][0]) / Lcell[i] );
    }
    
    //TODO: make for n dimensions
    if(numb_jk.size()==3){            
        int ID =
        jk_coords[2] * numb_jk[1] *numb_jk[2] +
        jk_coords[1] * numb_jk[1] +
        jk_coords[0];

        return ID;
    
    }else{
        std::cerr<< "#### ERROR [catalogue::jk_samples]: number of dimensions must be 3" << std::endl;
        exit (EXIT_FAILURE);
    }
    
}



// ==========================================================
// put objects into samples
// - samples build 3d grid mesh, spannig the the volume covered by the input catalogue
// - should be used for data in box. Light cone samplin still needs to be added
// ==========================================================
void catalogue::make_samples(std::vector < int > & numb_jk){

    
    //total number of samples
    int Nsamp=1;
    for(int i = 0; i < numb_jk.size(); i++){ Nsamp *= numb_jk[i]; }

    
    //initialize samples
    for(int i = 0; i < Nsamp; i++){
        sample s;
        samp.push_back(s);
    }
    
    
    //get limits of input cataloue
    std::vector < std::vector < double > > pos_limits = limits(input);
    //and add some tollerance to maxima to avoid that galaxies with
    //max position beeing put in JK sample with NK+1, which leads to segfault..
    for(int i=0; i<pos_limits.size(); i++){ pos_limits[i][1] *=1.00000001; }
    
    
    //length of sample cells along each dimension
    std::vector < double > Lcell;
    for(int i=0; i<pos_limits.size(); i++){
        Lcell.push_back( (pos_limits[i][1] - pos_limits[i][0]) / numb_jk[i] );
    }
    
    
    //add objects to samples
    for(int i=0; i<input.obj.size();i++){
        int ID = pos_to_ID(numb_jk, input.obj[i].pos, pos_limits, Lcell);
        samp[ID].obj.push_back(input.obj[i]);
    }
    
    
    //make sample center and edges    
    for(int i = 0; i < numb_jk[0]; i++){
        for(int j = 0; j < numb_jk[1]; j++){    
            for(int k = 0; k < numb_jk[2]; k++){
                
                std::vector <double> center = {
                    (i+0.5)*Lcell[0]-pos_limits[0][0],
                    (j+0.5)*Lcell[1]-pos_limits[1][0],
                    (k+0.5)*Lcell[2]-pos_limits[2][0]};
                
                int ID = pos_to_ID(numb_jk, center, pos_limits, Lcell);
                
                samp[ID].cent = center;
                
                //loop over edges
                for(int ie = 0; ie < 2; ie++){
                    for(int je = 0; je < 2; je++){
                        for(int ke = 0; ke < 2; ke++){
                                    
                            std::vector <double> coordinates = {
                                center[0] + (ie-0.5)*Lcell[0],
                                center[1] + (je-0.5)*Lcell[1],
                                center[2] + (ke-0.5)*Lcell[2]};
                            
                            samp[ID].edge.push_back(coordinates);

                        }
                    }
                }

                
            }
        }
    }
}



// ==========================================================
// free memory of input catalogue (use after samples are build)
// ==========================================================
void catalogue::delete_input(){
    
    input.obj.clear();
    input.cent.clear();
    input.edge.clear();

    sample().obj.swap(input.obj);
    sample().cent.swap(input.cent);
    sample().edge.swap(input.edge);
}



// ==========================================================
// generate random vector on a 3d spehere
// ==========================================================
std::vector<double> catalogue::rand_vec_sphere(double radius){

    
    int prcn=1000000;
    
    double r=1;
    double theta_rand = acos(((rand() % prcn)/double(prcn))*2-1);
    double phi_rand = ((rand() % prcn)/double(prcn))*2*M_PI;
    
    std::vector< double > vec;
    
    vec.push_back(radius * cos(phi_rand) * sin(theta_rand) );
    vec.push_back(radius * sin(phi_rand) * sin(theta_rand) );
    vec.push_back(radius * cos(theta_rand) );
    
    return vec;
    
}



// ==========================================================
// generate random vector on a 3d box
// ==========================================================
std::vector<double> catalogue::rand_vec_box(std::vector<double> x_lim, std::vector<double> y_lim, std::vector<double> z_lim){

    int prcn=1000000;
    
    double x = x_lim[0] + (rand() % prcn)/double(prcn)*(x_lim[1] - x_lim[0]);
    double y = y_lim[0] + (rand() % prcn)/double(prcn)*(y_lim[1] - y_lim[0]);
    double z = z_lim[0] + (rand() % prcn)/double(prcn)*(z_lim[1] - z_lim[0]);
    
    std::vector<double> vec = {x,y,z};
    
    return vec;
}



// ==========================================================
// generate random catalogue
// - objects with random positions and unit vectors with random positions
// ==========================================================
void catalogue::make_random(const parameters p){

    
    int prcn=10000;
    
    srand (p.rand_seed);
    
    for(int i=0; i<p.numb_rand; i++){
        
        double radius = 1.0;
        
        object obj_rand;
        
        obj_rand.pos = rand_vec_box(p.xlim_rand, p.ylim_rand, p.zlim_rand);
        obj_rand.vec_a = rand_vec_sphere(radius);
        obj_rand.vec_b = rand_vec_sphere(radius);
        
        random.obj.push_back(obj_rand);
    }
}



// ==========================================================
// write random catalogue
// - objects with random positions and unit vectors with random positions
// ==========================================================
void catalogue::write_input(const parameters p, const std::string filename){

    if(p.verbose>0)std::cout<<"# write: "<< filename <<std::endl;

    std::string d = p.delim_in;
    
    std::ofstream outdata;
    outdata.open(filename.c_str());
    outdata.setf(std::ios::fixed);
    outdata.precision (8);
    
    for(int i=0; i<input.obj.size(); i++){
        
        outdata <<
        input.obj[i].pos[0] <<d<<
        input.obj[i].pos[1] <<d<<
        input.obj[i].pos[2] <<d<<
        input.obj[i].vec_a[0] <<d<<
        input.obj[i].vec_a[1] <<d<<
        input.obj[i].vec_a[2] <<d<<
        input.obj[i].vec_b[0] <<d<<
        input.obj[i].vec_b[1] <<d<<
        input.obj[i].vec_b[2];
        if(i < input.obj.size()-1){ outdata<<std::endl; }
        
    }
    
}
