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

int main(int argc, char **argv)
{

    //===============================
    // read parameters
    //===============================

    if (argc < 2 || argc > 5)
    {
        std::cerr << "# USAGE: ./covo [parameter file] [option 1: catalogue_1] [option 2: catalogue_2]" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string fname_param = argv[1];

    parameters p;

    p.read(fname_param);

    if (argc >= 3)
    {
        p.fname_cat_1 = argv[2];
    }

    if (argc == 3)
    {
        p.fname_cat_2 = p.fname_cat_1;
    }

    if (argc >= 4)
    {
        p.fname_cat_2 = argv[3];
    }

    if (argc >= 5)
    {
        p.fname_out_suffix = argv[4];
    }

    if (p.verbose > 1)
        p.print();

    if (!p.check())
        exit(EXIT_FAILURE);

    //===============================
    // set up input catalogues
    //===============================

    catalogue cat_1;
    catalogue cat_2;

    if (p.make_rand)
    {
        if (p.verbose > 1)
            std::cout << std::endl
                      << "# ==== make random catalogue ====" << std::endl;

        if (p.mode == "box")
            cat_1.make_random_box(p);
        if (p.mode == "shell")
            cat_1.make_random_shell(p);

        cat_1.input = cat_1.random;
        cat_1.write_input(p, p.fname_rand);
        cat_2.input = cat_1.input;

        if (p.verbose > 1)
        {
            std::cout << "# size of catalogue 1: " << cat_1.input.obj.size() << std::endl;
            std::cout << "# catalogue 2 = catalogue 1 " << std::endl;
        }
    }
    else
    {

        if (p.verbose > 1)
            std::cout << std::endl
                      << "# ==== read data ====" << std::endl;

        cat_1.read(p, p.fname_cat_1);

        if (p.verbose > 1)
            std::cout << "# size of catalogue 1: " << cat_1.input.obj.size() << std::endl;

        if (p.fname_cat_1 == p.fname_cat_2)
        {
            cat_2 = cat_1;
            if (p.verbose > 1)
                std::cout << "\n# catalogue 2 = catalogue 1 " << std::endl;
        }
        else
        {
            cat_2.read(p, p.fname_cat_2);
            if (p.verbose > 1)
                std::cout << "\n# size of catalogue 2: " << cat_2.input.obj.size() << std::endl;
        }
    }

    // find min/max positions
    cat_1.get_pos_limits(p);
    cat_2.get_pos_limits(p);

    if (p.verbose > 0 && p.auto_limits && !p.make_rand)
    {
        std::cout << "\n# min. / max. coordinates in catalogue 1" << std::endl;
        cat_1.show_pos_limits(p);
        std::cout << std::endl;

        if (p.fname_cat_1 != p.fname_cat_2)
        {
            std::cout << "# min. / max. coordinates in catalogue 2" << std::endl;
            cat_2.show_pos_limits(p);
            std::cout << std::endl;
        }
    }

    // cut out catalogs at limits defined in parameter file
    if (p.mode == "box" && !p.auto_limits)
    {
        cat_1.cut_input(p);
        cat_2.cut_input(p);
    }

    // check differences in volume limits between catalogue 1 and 2
    double diff_toll = 0.01; // tollerance in relative difference between box limits [%]
    bool lim_diff = vol_lim_diff(p, cat_1, cat_2, diff_toll);
    if (lim_diff)
    {
        if (p.verbose > 0)
            std::cout << "#### WARINING: >" << diff_toll * 100 << " % difference in box limits of catalog 1 and catalog 2 ####" << std::endl;

        if (p.periodic_box && p.mode == "box" && p.auto_limits)
        {
            if (p.verbose > 0)
                std::cout << " -> set periodic_box = false" << std::endl;
            p.periodic_box = false;
        }
    }

    // cut out overlap region between catalogs
    if (p.auto_limits && !p.make_rand && p.fname_cat_1 != p.fname_cat_2)
    {
        if (p.verbose > 1)
            std::cout << "# cut out overlap region between cat 1 and cat 2\n"
                      << std::endl;
        cut_overlap(p, cat_1, cat_2);
    }

    // set box limit parameters to actual values determined above
    if (p.mode == "box" && !p.make_rand)
    {
        p.x_lim = cat_1.pos_limits_cart[0];
        p.y_lim = cat_1.pos_limits_cart[1];
        p.z_lim = cat_1.pos_limits_cart[2];
    }

    // normalize input vectors
    cat_1.normalize_vectors();
    cat_2.normalize_vectors();

    // make jack-knife samples
    if (p.verbose > 1)
        std::cout << "# make subsamples" << std::endl;

    if (p.mode == "shell")
    {
        cat_1.make_samples_shell(p);
        cat_2.make_samples_shell(p);
    }
    if (p.mode == "box")
    {
        cat_1.make_samples_box(p);
        cat_2.make_samples_box(p);
    }
    if (cat_1.samp.size() != cat_2.samp.size())
    {
        std::cerr << "#### ERROR: different subsamples in cat 1 and cat 2 ####" << std::endl;
        exit(EXIT_FAILURE);
    }

    // delete original input catalogue to free memory
    cat_1.delete_input();
    cat_2.delete_input();

    if (p.verbose > 1)
    {
        std::cout << "# number of subsamples: " << cat_1.samp.size() << std::endl;

        std::cout << "# total number of objects in sub samples " << std::endl;
        std::cout << "# cat1: " << cat_1.numb_objects() << ", cat2: " << cat_2.numb_objects() << std::endl;

        std::cout << "# mean number of objects per sub sample " << std::endl;
        std::cout << "# cat1: " << double(cat_1.numb_objects()) / double(cat_1.samp.size())
                  << ", cat2: " << double(cat_2.numb_objects()) / double(cat_2.samp.size()) << "\n"
                  << std::endl;
    }

    //===============================
    if (p.verbose > 1)
        std::cout << std::endl
                  << "# ==== compute correlation ==== " << std::endl;
    //===============================
    correlation corr;

    corr.compute(p, cat_1, cat_2);

    // check if ouput directory exists, make it if it doensn't
    if (!directory_exists(p.verbose, p.dir_out))
        make_directory(p.verbose, p.dir_out);

    // make filename
    if (p.fname_out == "auto")
    {
        p.fname_out = p.dir_out + p.fname_out_prefix + "_" + p.fname_out_suffix + "." + p.fname_out_extention;
    }
    if (p.fname_out == "input")
    {
        p.fname_out = p.fname_out_suffix;
    }

    // write results
    corr.write(p, p.fname_out);

    return 0;
}
