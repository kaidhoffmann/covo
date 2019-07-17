// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#include <string>
#include <iostream>
#include <sys/stat.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>


//===========================================================
// check if directory exists
//===========================================================
//from https://stackoverflow.com/questions/4980815/c-determining-if-directory-not-a-file-exists-in-linux
// note: function rquires <string.h> or <cstring>, not <string>
bool directory_exists(const short verbose, const std::string pathname)
{
    bool dir_ex = false;
    
    struct stat info;
    
    if( stat( pathname.c_str(), &info ) != 0 ){
        if(verbose>2) std::cerr<<"##### [io_data: directory_exists] WARNING: can't find "<< pathname <<" #####" << std::endl;
    }else if( info.st_mode & S_IFDIR ){
        dir_ex=true;
    }else{
        if(verbose>2) std::cerr<<"##### [io_data: directory_exists] WARNING: "<< pathname << "is no directory #####" << std::endl;
    }
    return dir_ex;

}



//===========================================================
// make directory
//===========================================================
//from https://codeyarns.com/2014/08/07/how-to-create-directory-using-c-on-linux/
void make_directory(const short verbose, const std::string pathname){
    
    if(verbose>2) std::cout<<"##### [io_data: make_directory]: make "<< pathname << " #####" << std::endl;
    
    
    const int dir_err = mkdir(pathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    if (dir_err == -1)
    {
        std::cerr<<"#### WARNING [io_data: make_directory]: Error creating directory! ####"<< std::endl;
        exit(EXIT_FAILURE);
    }
}


