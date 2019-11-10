// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#ifndef CORRELATION_H
#define CORRELATION_H

#include <vector>
#include "parameters.h"
#include "catalogue.h"

class correlation{
    
    struct vars{
        std::vector < unsigned long > counts;
        
        std::vector < double > r12_v1a;
        std::vector < double > r12_v1b;
        std::vector < double > r12_v2a;
        std::vector < double > r12_v2b;
        
        std::vector < double > v1a_v2a;
        std::vector < double > v1b_v2b;
        std::vector < double > v1a_v2b;
        std::vector < double > v1b_v2a;
    };

    struct results{
        std::vector < std::vector < vars > > sums_samps;
        std::vector < vars > jack_knife;
        vars total;
        vars error;
    };
    
    vars sums_pairs(
        const parameters p,
        const std::vector < catalogue::object > & obj_1,
        const std::vector < catalogue::object > & obj_2);

    
    double min_samp_dist(
        const parameters p,
        const std::vector < std::vector < double > > & edge_1,
        const std::vector < std::vector < double > > & edge_2);
    
    
    void sums_for_sample_combinations(const parameters p, catalogue & cat_1, catalogue & cat_2);
    
    void sum_up_sample_combinations(const parameters p);
    
    void jk_sampling(const parameters p);

    void normalize_total(const parameters p);
    
    void normalize_jk(const parameters p);

    void jk_errors(const parameters p);
    
    void print_cell_num(const int DIM, const int i);

    
    public:
    
    std::vector < std::vector < vars > > sums_samps;

    results res;

    void compute(const parameters p, catalogue & cat_1, catalogue & cat_2);
    
    void write(const parameters p, const std::string filename);
    
};


#endif // CORRELATION_H
