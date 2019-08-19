//see
//https://healpix.jpl.nasa.gov/html/cxxsupport/cxxutils_8h-source.html
//https://healpix.jpl.nasa.gov/html/Healpix_cxx/healpix__base_8cc-source.html
//https://healpix.jpl.nasa.gov/html/Healpix_cxx/healpix__base_8h-source.html

#include <iostream>
#include <cmath>
#include <climits>
#include "healpix.h"
#include "toolbox.h"

short healpix::ctab[];
short healpix::utab[];


// ==========================================================
//
// ==========================================================
healpix::Tablefiller::Tablefiller(){
    for (int m=0; m<0x100; ++m)
    {
    ctab[m] =
        (m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
        | ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4);
    utab[m] =
        (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
        | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }
}


healpix::Tablefiller healpix::Filler;


const int healpix::jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
const int healpix::jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };



// ==========================================================
//
// ==========================================================
inline double healpix::fmodulo (double v1, double v2) const{
    using namespace std;
    if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
    double tmp=fmod(v1,v2)+v2;
    return (tmp==v2) ? 0. : tmp;
}



// ==========================================================
// Adjusts the object to nside
// ==========================================================
void healpix::SetNside (int nside){
    order_  = nside2order(nside);
    if(order_<0) std::cout<< "#### WARNING: SetNside: nside must be power of 2 for nested maps ####" <<std::endl;
    nside_  = nside;
    npface_ = nside_*nside_;
    ncap_   = (npface_-nside_)<<1;    
    npix_   = 12*npface_;
    fact2_  = 4./npix_;
    fact1_  = (nside_<<1)*fact2_;
}


// ==========================================================
//
// ==========================================================
int healpix::nside2order (const int nside) const{
  if(nside<1) std::cout<< "invalid value for Nside" <<std::endl;
  if ((nside)&(nside-1)) return -1;
  return ilog2(nside);
}


// ==========================================================
//
// ==========================================================
int healpix::xyf2nest (int ix, int iy, int face_num) const{
  return (face_num<<(2*order_)) +
      (utab[ix&0xff] | (utab[ix>>8]<<16)
    | (utab[iy&0xff]<<1) | (utab[iy>>8]<<17));
}


// ==========================================================
//
// ==========================================================
void healpix::nest2xyf(int pix, int &ix, int &iy, int &face_num) const{
    
    face_num = pix>>(2*order_);
    pix &= (npface_-1);
    int raw = (pix&0x5555) | ((pix&0x55550000)>>15);
    ix = ctab[raw&0xff] | (ctab[raw>>8]<<4);
    pix >>= 1;
    raw = (pix&0x5555) | ((pix&0x55550000)>>15);
    iy = ctab[raw&0xff] | (ctab[raw>>8]<<4);
}



// ==========================================================
//
// ==========================================================
int healpix::ang2pix_z_phi_nest (double z, double phi) const {
    
    double za = fabs(z);
    double tt = fmodulo(phi,twopi) * inv_halfpi; // in [0,4)

    int face_num, ix, iy;
    if (za<=twothird) // Equatorial region
    {

        double temp1 = nside_*(0.5+tt);
        double temp2 = nside_*(z*0.75);
        int jp = int(temp1-temp2); // index of  ascending edge line
        int jm = int(temp1+temp2); // index of descending edge line
        int ifp = jp >> order_;  // in {0,4}
        int ifm = jm >> order_;
        if (ifp == ifm)           // faces 4 to 7
            face_num = (ifp==4) ? 4: ifp+4;
        else if (ifp < ifm)       // (half-)faces 0 to 3
            face_num = ifp;
        else                      // (half-)faces 8 to 11
            face_num = ifm + 8;

        ix = jm & (nside_-1);
        iy = nside_ - (jp & (nside_-1)) - 1;
    }
    else // polar region, za > 2/3
        {
            
        int ntt = int(tt);
        if (ntt>=4) ntt=3;
        double tp = tt-ntt;
        double tmp = nside_*sqrt(3*(1-za));

        int jp = int(tp*tmp); // increasing edge line index
        int jm = int((1.0-tp)*tmp); // decreasing edge line index
        if (jp>=nside_) jp = nside_-1; // for points too close to the boundary
        if (jm>=nside_) jm = nside_-1;
        if (z >= 0)
            {
            face_num = ntt;  // in {0,3}
            ix = nside_ - jm - 1;
            iy = nside_ - jp - 1;
            }
        else
            {
            face_num = ntt + 8; // in {8,11}
            ix =  jp;
            iy =  jm;
            }
        }

    return xyf2nest(ix,iy,face_num);
}




// ==========================================================
// Returns the number of the pixel which contains the angular coordinates (using nested scheme)
// ==========================================================
int healpix::ang2pix_nest(double theta, double phi) const{
    return ang2pix_z_phi_nest (cos(theta), phi);
}




// ==========================================================
//
// ==========================================================
void healpix::pix2ang_z_phi_nest (int pix, double &z, double &phi) {
    
    
    int nl4 = nside_*4;

    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);

    int jr = (jrll[face_num]<<order_) - ix - iy - 1;

    int nr, kshift;
    if (jr<nside_)
      {
      nr = jr;
      z = 1 - nr*nr*fact2_;
      kshift = 0;
      }
    else if (jr > 3*nside_)
      {
      nr = nl4-jr;
      z = nr*nr*fact2_ - 1;
      kshift = 0;
      }
    else
      {
      nr = nside_;
      z = (2*nside_-jr)*fact1_;
      kshift = (jr-nside_)&1;
      }

    int jp = (jpll[face_num]*nr + ix -iy + 1 + kshift) / 2;
    if (jp>nl4) jp-=nl4;
    if (jp<1) jp+=nl4;

    phi = (jp-(kshift+1)*0.5)*(halfpi/nr);
    
}

// ==========================================================
//Returns the angular coordinates of the center of the pixel with number pix
// ==========================================================

void healpix::pix2ang_nest (int pix, double &theta, double &phi){
    double z;
    pix2ang_z_phi_nest (pix,z,phi);
    theta = acos(z);
}


// ==========================================================
//
// ==========================================================
void healpix::make_mask(int nside, const std::vector <double> theta_lim, const std::vector <double> phi_lim){
    
    SetNside(nside);
    

    // loop over pixels
    for(int pix=0; pix<npix_; pix++){

        cell c;

        c.ID = pix;
        
        // --- step one ---
        
        // get angle from pixel ID
        double theta, phi;    
        pix2ang_nest (pix, theta, phi);
        
        // check if angles are outside of limits
        if(theta_lim[0] <= theta && theta <= theta_lim[1] && phi_lim[0] <= phi && phi <= phi_lim[1]){
            c.masked = false;
        }else{
            c.masked = true;
        }
        
        
        // --- step two ---
        //test of healpixels lie partly outside of range by
        //looping over points on 4 border lines and mask the cells they hit
        
        int Nbin_ang = 10000;

        double dtheta = (theta_lim[1]-theta_lim[0]) / double(Nbin_ang);
        double dphi = (phi_lim[1]-phi_lim[0]) / double(Nbin_ang);

        double toll_theta = (theta_lim[1]-theta_lim[0]) / 1000.;
        double toll_phi = (phi_lim[1]-phi_lim[0]) / 1000.;

        
        // 1
        for(int i=0; i<Nbin_ang; i++){
            phi = phi_lim[0] - dphi;
            theta = theta_lim[0] + i*dtheta;
            int ID = ang2pix_nest(theta, phi);
            if(ID==c.ID){c.masked = true;
            }
        }
        
        
        // 2
        for(int i=0; i<Nbin_ang; i++){
            phi = phi_lim[1] + dphi;
            theta = theta_lim[0] + i*dtheta;
            int ID = ang2pix_nest(theta, phi);
            if(ID==c.ID){c.masked = true;}
        }
        
        
        // 3
        for(int i=0; i<Nbin_ang; i++){
            phi = phi_lim[0] + i*dphi;
            theta = theta_lim[0] - dtheta;
            int ID = ang2pix_nest(theta, phi);
            if(ID==c.ID){c.masked = true;}
        }
        
        
        // 4
        for(int i=0; i<Nbin_ang; i++){
            phi = phi_lim[0] + i*dphi;
            theta = theta_lim[1] + dtheta;
            int ID = ang2pix_nest(theta, phi);
            if(ID==c.ID){c.masked = true;}
        }

        cells.push_back(c);

    }

}
