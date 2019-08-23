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
#include "toolbox.h"
#include "healpix.h"



// ==========================================================
// normalize vectors in input catalogue
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
// find minima and maxima of positions in cartesian cooredinates
// ==========================================================
std::vector < std::vector < double > > catalogue::find_limits_cart(sample smp){

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
}



// ==========================================================
// find minima and maxima of positions in spherical cooredinates
// ==========================================================
std::vector < std::vector < double > > catalogue::find_limits_sphere(sample smp){

    std::vector < std::vector < double > > lim;

    //initilize min/max position with first position vector
    std::vector < double > init = cart_to_sphere(smp.obj[0].pos);
    
    //std::cout<<init[0]<<"\t"<<init[1]<<"\t"<<init[2]<<std::endl;
    
    for(int i=0; i< init.size(); i++){
        lim.push_back({init[i],init[i]});
    }

    //std::cout<<lim[2][0]<<"\t"<<lim[2][1]<<std::endl;
    
    //search for min/max
    for(int i=0; i< smp.obj.size(); i++){
        
        std::vector <double> pos_sphere = cart_to_sphere(smp.obj[i].pos);
        
        for(int j=0; j<pos_sphere.size(); j++){
            
            
            if(pos_sphere[j] < lim[j][0]){
                lim[j][0] = pos_sphere[j];
            }

            if(pos_sphere[j] > lim[j][1]){
                lim[j][1] = pos_sphere[j];
            }
        }  
    }
    
    return lim;
}



// ==========================================================
// obtain cell ID for object from its position for cartesian subsamples
// ==========================================================
int catalogue::pos_to_ID_cart(
    const std::vector < int > & numb_jk_cart,
    const std::vector < double > & pos,
    const std::vector < std::vector < double > > & pos_limits,
    const std::vector < double > & Lcell){

    int jk_coords[pos.size()]={0};

    for(int i=0; i< pos.size(); i++){
        
        jk_coords[i] = int( (pos[i] - pos_limits[i][0]) / Lcell[i] );
    }
    
    //TODO: make for n dimensions
    if(numb_jk_cart.size()==3){            
        int ID =
        jk_coords[2] * numb_jk_cart[1] *numb_jk_cart[2] +
        jk_coords[1] * numb_jk_cart[1] +
        jk_coords[0];

        return ID;
    
    }else{
        std::cerr<< "#### ERROR [catalogue::jk_samples]: number of dimensions must be 3" << std::endl;
        exit (EXIT_FAILURE);
    }
    
}



// ==========================================================
// find mit min/max positions of sample
// ==========================================================
void catalogue::get_pos_limits(const parameters p){
    if(p.auto_limits && !p.make_rand){
        if(p.mode == "box") pos_limits_cart = find_limits_cart(input);
        if(p.mode == "shell") pos_limits_sphere = find_limits_sphere(input);
    }else{
        pos_limits_cart = {p.x_lim, p.y_lim, p.z_lim};
        pos_limits_sphere = {p.r_lim, p.theta_lim, p.phi_lim};
    }
};



// ==========================================================
// show mit min/max positions of sample
// ==========================================================
void catalogue::show_pos_limits(const parameters p){

    double f = 180/M_PI;
    
    if(p.mode == "box"){
        std::cout <<"# "<< pos_limits_cart[0][0] <<" < x < "<< pos_limits_cart[0][1] <<std::endl;
        std::cout <<"# "<< pos_limits_cart[1][0] <<" < y < "<< pos_limits_cart[1][1] <<std::endl;
        std::cout <<"# "<< pos_limits_cart[2][0] <<" < z < "<< pos_limits_cart[2][1] <<std::endl;
    }
    
    if(p.mode == "shell"){
        std::cout <<"# "<< pos_limits_sphere[0][0] <<" < r < "<< pos_limits_sphere[0][1] <<std::endl;
        std::cout <<"# "<< pos_limits_sphere[1][0]*f <<" < theta < "<< pos_limits_sphere[1][1]*f <<std::endl;
        std::cout <<"# "<< pos_limits_sphere[2][0]*f <<" < phi < "<< pos_limits_sphere[2][1]*f <<std::endl;
    }
};



// ==========================================================
// cut overla region between two input catalogues
// ==========================================================
void cut_overlap(const parameters p, catalogue & cat_1, catalogue & cat_2){
    
    //get common minimi and maxima psitions in cat_1 and  cat_2
    for(int i=0; i<3; i++){
    
        if(p.mode == "box"){

            //minima in cartesian coordinates
            if( cat_1.pos_limits_cart[i][0] < cat_2.pos_limits_cart[i][0] ){
                cat_1.pos_limits_cart[i][0] = cat_2.pos_limits_cart[i][0];
            }else{
                cat_2.pos_limits_cart[i][0] = cat_1.pos_limits_cart[i][0];
            }

            //minima in cartesian coordinates
            if( cat_1.pos_limits_cart[i][1] > cat_2.pos_limits_cart[i][1] ){
                cat_1.pos_limits_cart[i][1] = cat_2.pos_limits_cart[i][1];
            }else{
                cat_2.pos_limits_cart[i][1] = cat_1.pos_limits_cart[i][1];
            }
        }
        
        if(p.mode == "shell"){

            //minima in spherical coordinates
            if( cat_1.pos_limits_sphere[i][0] < cat_2.pos_limits_sphere[i][0] ){
                cat_1.pos_limits_sphere[i][0] = cat_2.pos_limits_sphere[i][0];
            }else{
                cat_2.pos_limits_sphere[i][0] = cat_1.pos_limits_sphere[i][0];
            }

            //spherical in spherical coordinates
            if( cat_1.pos_limits_sphere[i][1] > cat_2.pos_limits_sphere[i][1] ){
                cat_1.pos_limits_sphere[i][1] = cat_2.pos_limits_sphere[i][1];
            }else{
                cat_2.pos_limits_sphere[i][1] = cat_1.pos_limits_sphere[i][1];
            }
        }
        
    }
    
    cat_1.cut_input(p);
    cat_2.cut_input(p);
};

// ==========================================================
// erase points outside of min/max positions 
// ==========================================================
void catalogue::cut_input(const parameters p){
    
    for(int i=0; i<input.obj.size();i++){
        
        bool erase_object = false;

        if(p.mode == "box"){
            
            std::vector < double > pos = input.obj[i].pos;

            for(int j=0; j < pos.size() && erase_object==false; j++){
                if(
                    (pos[j] <= pos_limits_cart[j][0]) ||
                    (pos[j] >= pos_limits_cart[j][1])
                    ){
                        erase_object = true;
                    }                
            }
        }

        if(p.mode == "shell"){
        
            std::vector < double > pos = cart_to_sphere(input.obj[i].pos);

            for(int j=0; j < pos.size() && erase_object==false; j++){
                if(
                    (pos[j] <= pos_limits_sphere[j][0]) ||
                    (pos[j] >= pos_limits_sphere[j][1])
                    ){
                        erase_object = true;
                    }                
            }
        }
        
        if(erase_object){
            input.obj.erase(input.obj.begin()+i);
            i--;
        }
        
    }
}


// ==========================================================
// put objects into samples
// - samples build using 3d grid mesh, spannig the the volume covered by the input catalogue
// - should be used for data in box.
// ==========================================================
void catalogue::make_samples_box(const parameters p){
    
    //total number of samples
    int Nsamp=1;
    for(int i = 0; i < p.numb_jk_cart.size(); i++){ Nsamp *= p.numb_jk_cart[i]; }

    
    //initialize samples
    for(int i = 0; i < Nsamp; i++){
        sample s;
        samp.push_back(s);
    }
    
    //limits of input cataloue in cartesian coordinates x, y, z
    std::vector < std::vector < double > > pos_limits = pos_limits_cart;
    
    
    //add some tollerance to maxima to avoid that galaxies with
    //max position beeing put in JK sample with NK+1, which leads to segfault..
    for(int i=0; i<pos_limits.size(); i++){ pos_limits[i][1] *=1.00000001; }
    
    
    //length of sample cells along each dimension
    std::vector < double > Lcell;
    for(int i=0; i<pos_limits.size(); i++){
        Lcell.push_back( (pos_limits[i][1] - pos_limits[i][0]) / p.numb_jk_cart[i] );
    }
    
    
    //add objects to samples
    for(int i=0; i<input.obj.size();i++){
        int ID = pos_to_ID_cart(p.numb_jk_cart, input.obj[i].pos, pos_limits, Lcell);
        samp[ID].obj.push_back(input.obj[i]);
    }
    
    
    //make sample center and edges    
    for(int i = 0; i < p.numb_jk_cart[0]; i++){
        for(int j = 0; j < p.numb_jk_cart[1]; j++){    
            for(int k = 0; k < p.numb_jk_cart[2]; k++){

                std::vector <double> center = {
                    (i+0.5)*Lcell[0] + pos_limits[0][0],
                    (j+0.5)*Lcell[1] + pos_limits[1][0],
                    (k+0.5)*Lcell[2] + pos_limits[2][0]};

                int ID = pos_to_ID_cart(p.numb_jk_cart, center, pos_limits, Lcell);
                
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
// put objects into samples
// - samples build using healpix mesh
// - should be used for data in sphere.
// ==========================================================
void catalogue::make_samples_shell(const parameters p){
    
    //limits of input cataloue in spherical coordinates r, theta, phi
    std::vector < std::vector < double > > pos_limits = pos_limits_sphere;
    
    //make healpix mask
    healpix hp;
    hp.make_mask(p.nside, pos_limits[1], pos_limits[2]);
    
    //add some tollerance to radial maxima to avoid that galaxies with
    //max position beeing put in JK sample with NK+1, which leads to segfault..
    pos_limits[0][1] *=1.00000001;

    //radial bin width
    double dr = (pos_limits[0][1] - pos_limits[0][0]) / double(p.nrad);

    //initialize samples, add radial bin and healpix ID    
    for(int i=0; i<hp.cells.size(); i++){//loop over healpix cells
        if(!hp.cells[i].masked){
            for(int j=0; j<p.nrad; j++){//loop over radial bins
                sample s;
                s.hp_ID = hp.cells[i].ID;
                s.bin_rad = j;
                samp.push_back(s);
            }
        }
    }
    
    
    //add objects to samples
    for(int i=0; i<input.obj.size();i++){
        
        std::vector < double > pos_sphere = cart_to_sphere(input.obj[i].pos);
       
        int bin_rad = int((pos_sphere[0] - pos_limits[0][0]) / dr );
    
        long hp_ID = hp.ang2pix_nest(pos_sphere[1], pos_sphere[2]);
     
        for(int samp_ID=0; samp_ID<samp.size(); samp_ID++){
            
            if(samp[samp_ID].bin_rad == bin_rad && samp[samp_ID].hp_ID==hp_ID){
                samp[samp_ID].obj.push_back(input.obj[i]);
            }
        }
    }
    
    //sample edges: min, max of xyz corrdinates (like drawing a box, which encloses the sample)
    for(int samp_ID=0; samp_ID<samp.size(); samp_ID++){
        
        if(samp[samp_ID].obj.size()>0){//loop over samples
        
            std::vector< std::vector <double> > pos_limits;
        
            //initilize
            for(int dim=0; dim<samp[samp_ID].obj[0].pos.size(); dim++){//loop over dimensions
                double pos_ini = samp[samp_ID].obj[0].pos[dim];
                pos_limits.push_back({pos_ini, pos_ini});
            }
        
            //find min and max position in each dimension
            for(int obj_ID=0; obj_ID < samp[samp_ID].obj.size(); obj_ID++){//loop over objects in sample
                            
                for(int dim=0; dim < pos_limits.size(); dim++){//loop over dimensions
                    
                    //min
                    if(samp[samp_ID].obj[obj_ID].pos[dim] < pos_limits[dim][0]){
                        pos_limits[dim][0] = samp[samp_ID].obj[obj_ID].pos[dim];
                    }

                    //max
                    if(samp[samp_ID].obj[obj_ID].pos[dim] > pos_limits[dim][1]){
                        pos_limits[dim][1] = samp[samp_ID].obj[obj_ID].pos[dim];
                    }
                }
            }
        
            //loop over edges
            for(int ie = 0; ie < 2; ie++){
                for(int je = 0; je < 2; je++){
                    for(int ke = 0; ke < 2; ke++){

                        std::vector <double> coordinates =
                        {pos_limits[0][ie], pos_limits[1][je], pos_limits[2][ke]};
                        
                        samp[samp_ID].edge.push_back(coordinates);

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
// generate random vector on a 3d spehere
// ==========================================================
std::vector<double> catalogue::rand_vec_shell(const std::vector<double> & intrsq_lim, const std::vector<double> & cos_theta_lim, const std::vector<double> & phi_lim){

    int prcn=1000000;
    
    //inverse transform sampling
    
    double intrsq = intrsq_lim[0] + (rand() % prcn)/double(prcn)*(intrsq_lim[1] - intrsq_lim[0]);
    double r = pow(3*intrsq,1/3.);
    
    double cos_theta_rand = cos_theta_lim[0] + (rand() % prcn)/double(prcn)*(cos_theta_lim[1] - cos_theta_lim[0]);
    double theta_rand = acos(cos_theta_rand);
    
    double phi_rand = phi_lim[0] + (rand() % prcn)/double(prcn)*(phi_lim[1] - phi_lim[0]);
    
    std::vector< double > vec;

    vec.push_back(r * cos(phi_rand) * sin(theta_rand) );
    vec.push_back(r * sin(phi_rand) * sin(theta_rand) );
    vec.push_back(r * cos(theta_rand) );
    
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
// generate random catalogue in box
// - objects with random cartesian positions and unit vectors with random directions
// ==========================================================
void catalogue::make_random_box(const parameters p){

    int prcn=100000;
    
    srand (p.rand_seed);
    
    for(int i=0; i<p.numb_rand; i++){
        
        double radius = 1.0;
        
        object obj_rand;
        
        obj_rand.pos = rand_vec_box(p.x_lim_rand, p.y_lim_rand, p.z_lim_rand);
        obj_rand.vec_a = rand_vec_sphere(radius);
        obj_rand.vec_b = rand_vec_sphere(radius);
        
        random.obj.push_back(obj_rand);
    }
}


// ==========================================================
// generate random catalogue in shell
// - objects with random cartesian positions and unit vectors with random directions
// ==========================================================
void catalogue::make_random_shell(const parameters p){

    int prcn=100000;
    
    double r = 1.0;
    
    srand (p.rand_seed);
    
    for(int i=0; i<p.numb_rand; i++){
        
        object obj_rand;
        
        std::vector <double> intrsq_lim = {pow(p.r_lim_rand[0],3)/3., pow(p.r_lim_rand[1],3)/3.};
       
        std::vector <double> cos_theta_lim = {cos(p.theta_lim_rand[0]), cos(p.theta_lim_rand[1])};
        
        //if(p.mode shell)
        obj_rand.pos = rand_vec_shell(intrsq_lim, cos_theta_lim, p.phi_lim_rand);
        obj_rand.vec_a = rand_vec_sphere(r);
        obj_rand.vec_b = rand_vec_sphere(r);
        
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



// ==========================================================
// get total number of objects in subsamples
// ==========================================================
int catalogue::numb_objects(){

    int N = 0;
    
    for(int i=0; i<samp.size(); i++){
        N += samp[i].obj.size();
    }
    
    return N;
}

