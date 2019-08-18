/*
 *  This file was copied and modified from healpix_base.h in Healpix_cxx.
 *  Copyright (C) 2003, 2004, 2005, 2006 Max-Planck-Society
 *  (see https://healpix.jpl.nasa.gov/html/Healpix_cxx/files.html)
 *  the original author is Martin Reinecke
*/


#ifndef HEALPIX_BASE_H
#define HEALPIX_BASE_H

#include <cmath>
#include <vector>


class healpix{
    
    const double pi=3.141592653589793238462643383279502884197;
    const double twopi=6.283185307179586476925286766559005768394;
    const double inv_twopi=1.0/twopi;
    const double fourpi=12.56637061435917295385057353311801153679;
    const double halfpi=1.570796326794896619231321691639751442099;
    const double inv_halfpi=0.6366197723675813430755350534900574;
    const double inv_sqrt4pi = 0.2820947917738781434740397257803862929220;

    const double ln2 = 0.6931471805599453094172321214581766;
    const double inv_ln2 = 1.4426950408889634073599246810018921;
    const double ln10 = 2.3025850929940456840179914546843642;

    const double onethird=1.0/3.0;
    const double twothird=2.0/3.0;
    const double fourthird=4.0/3.0;

    const double degr2rad=pi/180.0;
    const double rad2degr=180.0/pi;
    
    static short ctab[0x100], utab[0x100];
    
    static const int jrll[];
    static const int jpll[];

    
    // The order of the map; -1 for nonhierarchical map.
    int order_;
    
    // The N_side parameter of the map; 0 if not allocated.
    int nside_;
    int npface_, ncap_, npix_;
    double fact1_, fact2_;

    //Adjusts the object to nside
    void SetNside (int nside);

    class Tablefiller{
      public:
        Tablefiller();
      };
      
      
    static Tablefiller Filler;
    
    //Returns the remainder of the division v1/v2.
    //The result is non-negative.
    //v1 can be positive or negative; v2 must be positive.
    inline double fmodulo (double v1, double v2) const;
    
    
    //Returns the remainder of the division v1/v2.
    //The result is non-negative.
    //v1 can be positive or negative; \a v2 must be positive.
    template<typename I> inline I imodulo (I v1, I v2)
    { I v=v1%v2; return (v>=0) ? v : v+v2; }
   
   
    //Returns the largest integer \a n that fulfills \a 2^n<=arg.
    template<typename I> inline int ilog2 (I arg) const
    {
        #ifdef __GNUC__
        if (arg==0) return 0;
        if (sizeof(I)==sizeof(int))
            return 8*sizeof(int)-1-__builtin_clz(arg);
        if (sizeof(I)==sizeof(long))
            return 8*sizeof(long)-1-__builtin_clzl(arg);
        if (sizeof(I)==sizeof(long long))
            return 8*sizeof(long long)-1-__builtin_clzll(arg);
        #endif
        int res=0;
        while (arg > 0xFFFF) { res+=16; arg>>=16; }
        if (arg > 0x00FF) { res|=8; arg>>=8; }
        if (arg > 0x000F) { res|=4; arg>>=4; }
        if (arg > 0x0003) { res|=2; arg>>=2; }
        if (arg > 0x0001) { res|=1; }
    
        return res;
    }
   
   
    struct cell {
        long ID;
        bool masked;
    };
    
    
    
    /*
    
    //! Returns the angular coordinates of the center of the pixel with number pix.
     pointing pix2ang (I pix) const
       {
       double z, phi, sth;
       bool have_sth;
       pix2loc (pix,z,phi,sth,have_sth);
       return have_sth ? pointing(atan2(sth,z),phi) : pointing(acos(z),phi);
    }
    
    */


    public:
    
    
    std::vector < cell > cells;
    
    void make_mask(int nside, const std::vector <double> theta_lim, const std::vector <double> phi_lim);
    
    int nside2order (int nside) const;
        
    // Returns the number of the pixel which contains the angular coordinates (using nested scheme)
    int ang2pix_nest (double theta, double phi) const;
    int ang2pix_z_phi_nest (double z, double phi) const;
        
    // Returns the angular coordinates of the center of the pixel with number pix
    void pix2ang_nest (int pix, double &theta, double &phi);
    void pix2ang_z_phi_nest (int pix, double &z, double &phi);
    
    int xyf2nest (int ix, int iy, int face_num) const;
    void nest2xyf(int pix, int &ix, int &iy, int &face_num) const;


};

#endif
