//from https://healpix.jpl.nasa.gov/html/Healpix_cxx/healpix__base2_8cc-source.html

#include <iostream>
#include <cmath>
#include <climits>
#include "healpix.h"

short healpix::utab[];

healpix::Tablefiller::Tablefiller()
  {
  for (int m=0; m<0x100; ++m)
    {
    utab[m] =
         (m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
      | ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }
  }

healpix::Tablefiller healpix::Filler;

inline double healpix::fmodulo (double v1, double v2){
    using namespace std;
    if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
    double tmp=fmod(v1,v2)+v2;
    return (tmp==v2) ? 0. : tmp;
}


int healpix::nside2order (int nside)
  {
  if(nside<1) std::cout<< "invalid value for Nside" <<std::endl;
  if ((nside)&(nside-1)) return -1;
  return ilog2(nside);
}


int healpix::xyf2nest (int ix, int iy, int face_num, int order) const{
  return (face_num<<(2*order)) +
      (utab[ix&0xff] | (utab[ix>>8]<<16)
    | (utab[iy&0xff]<<1) | (utab[iy>>8]<<17));
}



int healpix::ang2pix_z_phi_nest (double z, double phi, int nside) {
//Nested
    
    double za = fabs(z);
    double tt = fmodulo(phi,twopi) * inv_halfpi; // in [0,4)

    int order = nside2order(nside);
    
    int npface = nside*nside;

    int ncap  = (npface-nside)<<1;
    int npix  = 12*npface;

    int face_num, ix, iy;
    if (za<=twothird) // Equatorial region
    {

        double temp1 = nside*(0.5+tt);
        double temp2 = nside*(z*0.75);
        int jp = int(temp1-temp2); // index of  ascending edge line
        int jm = int(temp1+temp2); // index of descending edge line
        int ifp = jp >> order;  // in {0,4}
        int ifm = jm >> order;
        if (ifp == ifm)           // faces 4 to 7
            face_num = (ifp==4) ? 4: ifp+4;
        else if (ifp < ifm)       // (half-)faces 0 to 3
            face_num = ifp;
        else                      // (half-)faces 8 to 11
            face_num = ifm + 8;

        ix = jm & (nside-1);
        iy = nside - (jp & (nside-1)) - 1;
    }
    else // polar region, za > 2/3
        {
            
        //std::cout<<"polar region "<<z<<"\t"<<za<<"\t"<<phi<<"\t"<<twothird<<std::endl;

        int ntt = int(tt);
        if (ntt>=4) ntt=3;
        double tp = tt-ntt;
        double tmp = nside*sqrt(3*(1-za));

        int jp = int(tp*tmp); // increasing edge line index
        int jm = int((1.0-tp)*tmp); // decreasing edge line index
        if (jp>=nside) jp = nside-1; // for points too close to the boundary
        if (jm>=nside) jm = nside-1;
        if (z >= 0)
            {
            face_num = ntt;  // in {0,3}
            ix = nside - jm - 1;
            iy = nside - jp - 1;
            }
        else
            {
            face_num = ntt + 8; // in {8,11}
            ix =  jp;
            iy =  jm;
            }
        }

    return xyf2nest(ix,iy,face_num, order);
}

    
int healpix::ang2pix_nest(double theta, double phi, int nside){
    return ang2pix_z_phi_nest (cos(theta), phi, nside);
}
