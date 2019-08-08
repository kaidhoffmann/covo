// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019


#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "directory.h"
#include "parameters.h"
#include "catalogue.h"
#include "correlation.h"


int main(int argc, char **argv){
    
    
    //===============================
    // read parameters
    //===============================
    
    if(argc < 2 || argc > 5){
        std::cerr<<"# USAGE: ./covo [parameter file] [option 1: catalogue_1] [option 2: catalogue_2]"<<std::endl;
        exit (EXIT_FAILURE);
    }
    
    std::string fname_param = argv[1];
    
    parameters p;
    
    p.read(fname_param);
            
    if(argc>=3){
        p.fname_cat_1 = argv[2];
    }
        
    if(argc==3){
        p.fname_cat_2 = p.fname_cat_1;
    }

    if(argc>=4){
        p.fname_cat_2 = argv[3];
    }
    
    if(argc>=5){
        p.fname_out_suffix = argv[4];        
    }
    
    if(p.verbose > 1) p.print();
    
    if(!p.check()) exit (EXIT_FAILURE);

    
    
    //===============================
    // set up input catalogues
    //===============================
    catalogue cat_1;
    catalogue cat_2;
    
    if(p.make_rand){
        if(p.verbose > 1) std::cout<<std::endl<<"# ==== make random catalogue ====" << std::endl;
        
        cat_1.make_random(p);
        cat_1.input = cat_1.random;
        cat_1.write_input(p, p.fname_rand);
        cat_2.input = cat_1.input;
        
        if(p.verbose > 1){
            std::cout<<"# size of catalogue 1: "<<cat_1.input.obj.size()<<std::endl;
            std::cout<<"# catalogue 2 = catalogue 1 "<<std::endl;
        }
        
    }else{
            
        
        if(p.verbose > 1) std::cout<<std::endl<<"# ==== read data ====" << std::endl;

        cat_1.read(p, p.fname_cat_1);
        
        if(p.verbose > 1) std::cout<<"# size of catalogue 1: "<<cat_1.input.obj.size()<<std::endl;
    
        if(p.fname_cat_1==p.fname_cat_2){
            cat_2 = cat_1;
            if(p.verbose > 1) std::cout<<"# catalogue 2 = catalogue 1 "<<std::endl;
        }else{
            cat_2.read(p, p.fname_cat_2);
            if(p.verbose > 1) std::cout<<"# size of catalogue 2: "<<cat_2.input.obj.size()<<std::endl;        
        }
    
    }
    //normalize input vectors
    cat_1.normalize_vectors();
    cat_2.normalize_vectors();
    
    //make jack-knife samples
    cat_1.make_samples_cart(p.numb_jk);
    cat_2.make_samples_cart(p.numb_jk);
    
    //healpix
    //cat_1.make_samples_healpix(4);
    
    
    //delete original input catalogue to free memory
    cat_1.delete_input();
    cat_2.delete_input();
    
    
    
    //===============================
    if(p.verbose > 1) std::cout<<std::endl<<"# ==== compute correlation ==== " << std::endl;
    //===============================
    correlation corr;

    corr.compute(p, cat_1, cat_2);
    
    //check if ouput directory exists, make it if it doensn't
    if(!directory_exists(p.verbose, p.dir_out)) make_directory(p.verbose, p.dir_out);
    
    //make filename
    if(p.fname_out == "auto"){
        p.fname_out = p.dir_out + p.fname_out_prefix + "_" + p.fname_out_suffix + "." + p.fname_out_extention;
    }
    
    //write results
    corr.write(p, p.fname_out);
    
    
    
    return 0;
    
}
