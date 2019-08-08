// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#ifndef TOOLBOX_H
#define TOOLBOX_H


#include <vector>
#include <cmath>



// absolute value of vector
double vabs(const std::vector<double> vec);


// convert cartesian to spehrical coordinates
std::vector<double> cart_to_sphere(const std::vector<double> pos_cart);


// convert spehrical coordinates to cartesian
std::vector<double> sphere_to_cart(const std::vector<double> pos_sphere);


#endif // TOOLBOX_H
