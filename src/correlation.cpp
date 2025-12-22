// Kai Hoffmann
// Institute for Computational Science (ICS)
// University of Zurich
// 2019

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <functional> // std::divides
#include <array>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "parameters.h"
#include "catalogue.h"
#include "correlation.h"
#include "toolbox.h"

namespace
{
inline double dot_arr_vec(const double *a, const std::vector<double> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double dot_vec_vec(const std::vector<double> &a, const std::vector<double> &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double pow_exp(double v, double expip)
{
    double av = std::abs(v);
    if (expip == 1.0)
        return av;
    if (expip == 2.0)
        return av * av;
    return pow(av, expip);
}

inline void print_progress_bar(int completed, int total)
{
    const int bar_width = 50;
    const double percent = 100.0 * completed / total;
    const int filled = static_cast<int>(bar_width * completed / total);
    
    std::cout << "\r# [";
    for (int i = 0; i < bar_width; i++)
    {
        std::cout << (i < filled ? '=' : ' ');
    }
    std::cout << "] " << static_cast<int>(percent) << "% (" << completed << "/" << total << ")";
    std::cout.flush();
}

template<typename CounterType>
inline void update_and_print_progress(CounterType& completed, int total)
{
    int count = ++completed;
#ifdef _OPENMP
#pragma omp critical
    {
        print_progress_bar(count, total);
        if (count == total)
            std::cout << std::endl;
    }
#else
    print_progress_bar(count, total);
    if (count == total)
        std::cout << std::endl;
#endif
}
} // namespace

// ==========================================================
// mimimal distance between two samples defined by their edges
// ==========================================================
double correlation::min_samp_dist(
    const parameters p,
    const std::vector<std::vector<double>> &edge_1,
    const std::vector<std::vector<double>> &edge_2)
{

    double r_max = float(p.r_max);
    bool periodic = false;

    if (p.mode == "box" && p.periodic_box)
        periodic = true;

    int dim = edge_1[0].size();

    double Lbox[dim] = {0};
    if (p.mode == "box")
    {
        Lbox[0] = p.x_lim[1] - p.x_lim[0];
        Lbox[1] = p.y_lim[1] - p.y_lim[0];
        Lbox[2] = p.z_lim[1] - p.z_lim[0];
    }

    // initialize min distance
    double dist_min = 0;
    for (int k = 0; k < dim; k++)
    {
        double dk = edge_2[0][k] - edge_1[0][k];
        if (periodic)
        {
            dk = periodic_distance(dk, Lbox[k], r_max);
        }
        dist_min += dk * dk;
    }
    dist_min = pow(dist_min, 0.5);

    // find min distances between edges of cell pair
    for (int i = 0; i < edge_1.size(); i++)
    {
        for (int j = 0; j < edge_2.size(); j++)
        {

            double dist = 0;
            for (int k = 0; k < dim; k++)
            {
                double dk = edge_2[j][k] - edge_1[i][k];
                if (periodic)
                {
                    dk = periodic_distance(dk, Lbox[k], r_max);
                }
                dist += dk * dk;
            }
            dist = pow(dist, 0.5);
            if (dist < dist_min)
            {
                dist_min = dist;
            }
        }
    }
    return dist_min;
};

// ==========================================================
// finds pairs of objects belonging to two catalogues
// and sum up counts and inner products in bins of pair distance
// ==========================================================
correlation::vars correlation::sums_pairs(
    const parameters p,
    const std::vector<catalogue::object> &obj_1,
    const std::vector<catalogue::object> &obj_2)
{

    const double dr = (p.r_max - p.r_min) / static_cast<double>(p.numb_bin);
    const double r_max = static_cast<double>(p.r_max);
    const double r_min_sq = p.r_min * p.r_min;
    const double r_max_sq = p.r_max * p.r_max;
    const double lg_r_min = log10(p.r_min);
    const double dlg_r = (log10(p.r_max) - log10(p.r_min)) / static_cast<double>(p.numb_bin);
    const double inv_dr = 1.0 / dr;
    const double inv_dlg_r = 1.0 / dlg_r;
    const double expip = p.expip;

    vars sums_samp;

    // initialize
    sums_samp.counts.assign(p.numb_bin, 0);
    sums_samp.r12_v1a.clear();
    sums_samp.r12_v1b.clear();
    sums_samp.r12_v2a.clear();
    sums_samp.r12_v2b.clear();

    sums_samp.v1a_v2a.clear();
    sums_samp.v1b_v2b.clear();
    sums_samp.v1a_v2b.clear();
    sums_samp.v1b_v2a.clear();

    if (p.r12_v1a)
    {
        sums_samp.r12_v1a.resize(p.numb_bin, 0.0);
    }
    if (p.r12_v1b)
    {
        sums_samp.r12_v1b.resize(p.numb_bin, 0.0);
    }
    if (p.r12_v2a)
    {
        sums_samp.r12_v2a.resize(p.numb_bin, 0.0);
    }
    if (p.r12_v2b)
    {
        sums_samp.r12_v2b.resize(p.numb_bin, 0.0);
    }

    if (p.v1a_v2a)
    {
        sums_samp.v1a_v2a.resize(p.numb_bin, 0.0);
    }
    if (p.v1b_v2b)
    {
        sums_samp.v1b_v2b.resize(p.numb_bin, 0.0);
    }
    if (p.v1a_v2b)
    {
        sums_samp.v1a_v2b.resize(p.numb_bin, 0.0);
    }
    if (p.v1b_v2a)
    {
        sums_samp.v1b_v2a.resize(p.numb_bin, 0.0);
    }

    // box size
    std::array<double, 3> Lbox{};
    if (p.mode == "box")
    {
        Lbox[0] = p.x_lim[1] - p.x_lim[0];
        Lbox[1] = p.y_lim[1] - p.y_lim[0];
        Lbox[2] = p.z_lim[1] - p.z_lim[0];
    }

    const bool periodic = (p.mode == "box" && p.periodic_box);

    const std::size_t n1 = obj_1.size();
    const std::size_t n2 = obj_2.size();

    for (std::size_t i = 0; i < n1; i++)
    {
        for (std::size_t j = 0; j < n2; j++)
        {

            std::array<double, 3> d{};

            d[0] = obj_2[j].pos[0] - obj_1[i].pos[0];
            if (periodic)
            {
                d[0] = periodic_distance(d[0], Lbox[0], r_max);
            }

            if (std::abs(d[0]) < r_max)
            {

                d[1] = obj_2[j].pos[1] - obj_1[i].pos[1];
                if (periodic)
                {
                    d[1] = periodic_distance(d[1], Lbox[1], r_max);
                }

                if (std::abs(d[1]) < r_max)
                {

                    d[2] = obj_2[j].pos[2] - obj_1[i].pos[2];
                    if (periodic)
                    {
                        d[2] = periodic_distance(d[2], Lbox[2], r_max);
                    }

                    if (std::abs(d[2]) < r_max)
                    {

                        double r_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

                        if (r_sq < r_max_sq)
                        {
                            if (r_sq > r_min_sq)
                            {

                                int bin;

                                if (p.lg_bins)
                                {
                                    bin = static_cast<int>((0.5 * log10(r_sq) - lg_r_min) * inv_dlg_r);
                                }
                                else
                                {
                                    const double r_abs = sqrt(r_sq);
                                    bin = static_cast<int>((r_abs - p.r_min) * inv_dr);
                                }

                                double r_inv = 1.0f / sqrt(r_sq);
                                double r_vec[3] = {d[0] * r_inv, d[1] * r_inv, d[2] * r_inv};

                                sums_samp.counts[bin]++;

                                if (p.r12_v1a)
                                {
                                    sums_samp.r12_v1a[bin] += pow_exp(dot_arr_vec(r_vec, obj_1[i].vec_a), expip);
                                }

                                if (p.r12_v1b)
                                {
                                    sums_samp.r12_v1b[bin] += pow_exp(dot_arr_vec(r_vec, obj_1[i].vec_b), expip);
                                }

                                if (p.r12_v2a)
                                {
                                    sums_samp.r12_v2a[bin] += pow_exp(dot_arr_vec(r_vec, obj_2[j].vec_a), expip);
                                }

                                if (p.r12_v2b)
                                {
                                    sums_samp.r12_v2b[bin] += pow_exp(dot_arr_vec(r_vec, obj_2[j].vec_b), expip);
                                }

                                if (p.v1a_v2a)
                                {
                                    sums_samp.v1a_v2a[bin] += pow_exp(dot_vec_vec(obj_1[i].vec_a, obj_2[j].vec_a), expip);
                                }

                                if (p.v1b_v2b)
                                {
                                    sums_samp.v1b_v2b[bin] += pow_exp(dot_vec_vec(obj_1[i].vec_b, obj_2[j].vec_b), expip);
                                }

                                if (p.v1a_v2b)
                                {
                                    sums_samp.v1a_v2b[bin] += pow_exp(dot_vec_vec(obj_1[i].vec_a, obj_2[j].vec_b), expip);
                                }

                                if (p.v1b_v2a)
                                {
                                    sums_samp.v1b_v2a[bin] += pow_exp(dot_vec_vec(obj_1[i].vec_b, obj_2[j].vec_a), expip);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return sums_samp;
}

// ==========================================================
// search for pairs in all samples with distance r < r_max
// ==========================================================
void correlation::sums_for_sample_combinations(const parameters p, catalogue &cat_1, catalogue &cat_2)
{
    const std::size_t n1 = cat_1.samp.size();
    const std::size_t n2 = cat_2.samp.size();
    sums_samps.resize(n1);

    const int total = static_cast<int>(n1);

#ifdef _OPENMP
    if (p.num_threads > 0)
    {
        omp_set_num_threads(p.num_threads);
    }
    std::atomic<int> completed(0);
#pragma omp parallel for schedule(dynamic)
#else
    int completed = 0;
#endif
    for (int i = 0; i < total; i++)
    {
        if (p.verbose > 0)
        {
            update_and_print_progress(completed, total);
        }

        std::vector<vars> sums_i;
        for (std::size_t j = 0; j < n2; j++)
        {

            double dist_samps = min_samp_dist(p, cat_1.samp[i].edge, cat_2.samp[j].edge);

            if (dist_samps <= p.r_max)
            {
                vars sums_ij = sums_pairs(p, cat_1.samp[i].obj, cat_2.samp[j].obj);
                sums_i.push_back(sums_ij);
            }
        }

        sums_samps[i] = sums_i;
    }
}

// ==========================================================
// sum up results for pairs in all sample combinations
// -gives results in total volume
// ==========================================================
void correlation::sum_up_sample_combinations(const parameters p)
{

    res.total.counts.clear();

    res.total.r12_v1a.clear();
    res.total.r12_v1b.clear();
    res.total.r12_v2a.clear();
    res.total.r12_v2b.clear();

    res.total.v1a_v2a.clear();
    res.total.v1b_v2b.clear();
    res.total.v1a_v2b.clear();
    res.total.v1b_v2a.clear();

    res.total.counts.resize(p.numb_bin, 0);

    res.total.r12_v1a.resize(p.numb_bin, 0.0);
    res.total.r12_v1b.resize(p.numb_bin, 0.0);
    res.total.r12_v2a.resize(p.numb_bin, 0.0);
    res.total.r12_v2b.resize(p.numb_bin, 0.0);

    res.total.v1a_v2a.resize(p.numb_bin, 0.0);
    res.total.v1b_v2b.resize(p.numb_bin, 0.0);
    res.total.v1a_v2b.resize(p.numb_bin, 0.0);
    res.total.v1b_v2a.resize(p.numb_bin, 0.0);

    for (int i = 0; i < sums_samps.size(); i++)
    {
        for (int j = 0; j < sums_samps[i].size(); j++)
        {

            sum_vectors(res.total.counts, sums_samps[i][j].counts);

            if (p.r12_v1a)
            {
                sum_vectors(res.total.r12_v1a, sums_samps[i][j].r12_v1a);
            }

            if (p.r12_v1b)
            {
                sum_vectors(res.total.r12_v1b, sums_samps[i][j].r12_v1b);
            }

            if (p.r12_v2a)
            {
                sum_vectors(res.total.r12_v2a, sums_samps[i][j].r12_v2a);
            }

            if (p.r12_v2b)
            {
                sum_vectors(res.total.r12_v2b, sums_samps[i][j].r12_v2b);
            }

            if (p.v1a_v2a)
            {
                sum_vectors(res.total.v1a_v2a, sums_samps[i][j].v1a_v2a);
            }

            if (p.v1b_v2b)
            {
                sum_vectors(res.total.v1b_v2b, sums_samps[i][j].v1b_v2b);
            }

            if (p.v1a_v2b)
            {
                sum_vectors(res.total.v1a_v2b, sums_samps[i][j].v1a_v2b);
            }

            if (p.v1b_v2a)
            {
                sum_vectors(res.total.v1b_v2a, sums_samps[i][j].v1b_v2a);
            }
        }
    }
}

// ==========================================================
// get results in "neglect one" jack-knife samples by substracting
// sample combinations involving a given jk cell from the total result
// ==========================================================
void correlation::jk_sampling(const parameters p)
{

    int numb_samp = sums_samps.size();

    res.jack_knife.resize(numb_samp, res.total);

    for (int s = 0; s < numb_samp; s++)
    {
        for (int i = 0; i < sums_samps.size(); i++)
        {
            for (int j = 0; j < sums_samps[i].size(); j++)
            {
                if (i == s || j == s)
                {

                    subtract_vectors(res.jack_knife[s].counts, sums_samps[i][j].counts);

                    if (p.r12_v1a)
                    {
                        subtract_vectors(res.jack_knife[s].r12_v1a, sums_samps[i][j].r12_v1a);
                    }

                    if (p.r12_v1b)
                    {
                        subtract_vectors(res.jack_knife[s].r12_v1b, sums_samps[i][j].r12_v1b);
                    }

                    if (p.r12_v2a)
                    {
                        subtract_vectors(res.jack_knife[s].r12_v2a, sums_samps[i][j].r12_v2a);
                    }

                    if (p.r12_v2b)
                    {
                        subtract_vectors(res.jack_knife[s].r12_v2b, sums_samps[i][j].r12_v2b);
                    }

                    if (p.v1a_v2a)
                    {
                        subtract_vectors(res.jack_knife[s].v1a_v2a, sums_samps[i][j].v1a_v2a);
                    }

                    if (p.v1b_v2b)
                    {
                        subtract_vectors(res.jack_knife[s].v1b_v2b, sums_samps[i][j].v1b_v2b);
                    }

                    if (p.v1a_v2b)
                    {
                        subtract_vectors(res.jack_knife[s].v1a_v2b, sums_samps[i][j].v1a_v2b);
                    }

                    if (p.v1b_v2a)
                    {
                        subtract_vectors(res.jack_knife[s].v1b_v2a, sums_samps[i][j].v1b_v2a);
                    }
                }
            }
        }
    }
}

// ==========================================================
// get average results for pairs in full volume by normalizin with counts
// ==========================================================
void correlation::normalize_total(const parameters p)
{

    // normalize total
    std::vector<double> norm(res.total.counts.begin(), res.total.counts.end());

    divide_vectors(res.total.r12_v1a, norm);
    divide_vectors(res.total.r12_v1b, norm);
    divide_vectors(res.total.r12_v2a, norm);
    divide_vectors(res.total.r12_v2b, norm);

    divide_vectors(res.total.v1a_v2a, norm);
    divide_vectors(res.total.v1b_v2b, norm);
    divide_vectors(res.total.v1a_v2b, norm);
    divide_vectors(res.total.v1b_v2a, norm);
}

// ==========================================================
// get average results for pairs in jk samples volume by normalizin with counts
// ==========================================================
void correlation::normalize_jk(const parameters p)
{

    // normalize jack-knife samples
    int numb_samp = sums_samps.size();

    for (int s = 0; s < numb_samp; s++)
    {

        std::vector<double> norm(res.jack_knife[s].counts.begin(), res.jack_knife[s].counts.end());

        divide_vectors(res.jack_knife[s].r12_v1a, norm);
        divide_vectors(res.jack_knife[s].r12_v1b, norm);
        divide_vectors(res.jack_knife[s].r12_v2a, norm);
        divide_vectors(res.jack_knife[s].r12_v2b, norm);

        divide_vectors(res.jack_knife[s].v1a_v2a, norm);
        divide_vectors(res.jack_knife[s].v1b_v2b, norm);
        divide_vectors(res.jack_knife[s].v1a_v2b, norm);
        divide_vectors(res.jack_knife[s].v1b_v2a, norm);
    }
}

// ==========================================================
// get jack-knife estimates for standard deviation
// ==========================================================
void correlation::jk_errors(const parameters p)
{

    int numb_samp = sums_samps.size();

    res.error.counts.clear();

    res.error.r12_v1a.clear();
    res.error.r12_v1b.clear();
    res.error.r12_v2a.clear();
    res.error.r12_v2b.clear();

    res.error.v1a_v2a.clear();
    res.error.v1b_v2b.clear();
    res.error.v1a_v2b.clear();
    res.error.v1b_v2a.clear();

    res.error.counts.resize(p.numb_bin, 0);

    res.error.r12_v1a.resize(p.numb_bin, 0.0);
    res.error.r12_v1b.resize(p.numb_bin, 0.0);
    res.error.r12_v2a.resize(p.numb_bin, 0.0);
    res.error.r12_v2b.resize(p.numb_bin, 0.0);

    res.error.v1a_v2a.resize(p.numb_bin, 0.0);
    res.error.v1b_v2b.resize(p.numb_bin, 0.0);
    res.error.v1a_v2b.resize(p.numb_bin, 0.0);
    res.error.v1b_v2a.resize(p.numb_bin, 0.0);

    for (int s = 0; s < numb_samp; s++)
    {
        for (int i = 0; i < p.numb_bin; i++)
        {

            res.error.counts[i] += pow(res.total.counts[i] - res.jack_knife[s].counts[i], 2) * double(numb_samp - 1) / double(numb_samp);

            if (p.r12_v1a)
            {
                res.error.r12_v1a[i] += pow(res.total.r12_v1a[i] - res.jack_knife[s].r12_v1a[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.r12_v1b)
            {
                res.error.r12_v1b[i] += pow(res.total.r12_v1b[i] - res.jack_knife[s].r12_v1b[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.r12_v2a)
            {
                res.error.r12_v2a[i] += pow(res.total.r12_v2a[i] - res.jack_knife[s].r12_v2a[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.r12_v2b)
            {
                res.error.r12_v2b[i] += pow(res.total.r12_v2b[i] - res.jack_knife[s].r12_v2b[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.v1a_v2a)
            {
                res.error.v1a_v2a[i] += pow(res.total.v1a_v2a[i] - res.jack_knife[s].v1a_v2a[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.v1b_v2b)
            {
                res.error.v1b_v2b[i] += pow(res.total.v1b_v2b[i] - res.jack_knife[s].v1b_v2b[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.v1a_v2b)
            {
                res.error.v1a_v2b[i] += pow(res.total.v1a_v2b[i] - res.jack_knife[s].v1a_v2b[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }

            if (p.v1b_v2a)
            {
                res.error.v1b_v2a[i] += pow(res.total.v1b_v2a[i] - res.jack_knife[s].v1b_v2a[i], 2) * double(numb_samp - 1) / double(numb_samp);
            }
        }
    }

    for (int i = 0; i < p.numb_bin; i++)
    {

        res.error.counts[i] = sqrt(res.error.counts[i]);

        if (p.r12_v1a)
        {
            res.error.r12_v1a[i] = sqrt(res.error.r12_v1a[i]);
        }

        if (p.r12_v1b)
        {
            res.error.r12_v1b[i] = sqrt(res.error.r12_v1b[i]);
        }

        if (p.r12_v2a)
        {
            res.error.r12_v2a[i] = sqrt(res.error.r12_v2a[i]);
        }

        if (p.r12_v2b)
        {
            res.error.r12_v2b[i] = sqrt(res.error.r12_v2b[i]);
        }

        if (p.v1a_v2a)
        {
            res.error.v1a_v2a[i] = sqrt(res.error.v1a_v2a[i]);
        }

        if (p.v1b_v2b)
        {
            res.error.v1b_v2b[i] = sqrt(res.error.v1b_v2b[i]);
        }

        if (p.v1a_v2b)
        {
            res.error.v1a_v2b[i] = sqrt(res.error.v1a_v2b[i]);
        }

        if (p.v1b_v2a)
        {
            res.error.v1b_v2a[i] = sqrt(res.error.v1b_v2a[i]);
        }
    }
}

// ==========================================================
// driver function for all the steps of computing the correlation
// ==========================================================
void correlation::compute(const parameters p, catalogue &cat_1, catalogue &cat_2)
{

    time_t start;
    time(&start);

    // get results for combinations of samples
    sums_for_sample_combinations(p, cat_1, cat_2);

    // get total sums
    sum_up_sample_combinations(p);

    // subtract results in jack-knife regions
    jk_sampling(p);

    // normalize total
    normalize_total(p);

    // normalize jk
    normalize_jk(p);

    // normalize jk
    jk_errors(p);

    time_t end;
    time(&end);

    if (p.verbose > 0)
    {
        std::cout << std::endl
                  << "# " << double(difftime(end, start)) << " seconds for computing correlation" << std::endl;
    }
}

// ==========================================================
// write final to output file
// ==========================================================
void correlation::write(const parameters p, const std::string filename)
{

    if (p.verbose > 0)
        std::cout << "# write: " << filename << std::endl;

    double dr = (p.r_max - p.r_min) / double(p.numb_bin);
    double dlg_r = (log10(p.r_max) - log10(p.r_min)) / double(p.numb_bin);

    std::string d = p.delim_out;

    std::ofstream outdata;
    outdata.open(filename.c_str());
    outdata.setf(std::ios::fixed);
    outdata.precision(8);

    if (p.header_out)
    {

        outdata << "# ";
        outdata << "bin";
        outdata << d << "r";
        outdata << d << "counts";

        if (p.r12_v1a)
            outdata << d << "r12_v1a" << d << "r12_v1a_std";
        if (p.r12_v1b)
            outdata << d << "r12_v1b" << d << "r12_v1b_std";
        if (p.r12_v2a)
            outdata << d << "r12_v2a" << d << "r12_v2a_std";
        if (p.r12_v2b)
            outdata << d << "r12_v2b" << d << "r12_v2b_std";

        if (p.v1a_v2a)
            outdata << d << "v1a_v2a" << d << "v1a_v2a_std";
        if (p.v1b_v2b)
            outdata << d << "v1b_v2b" << d << "v1b_v2b_std";
        if (p.v1a_v2b)
            outdata << d << "v1a_v2b" << d << "v1a_v2b_std";
        if (p.v1b_v2a)
            outdata << d << "v1b_v2a" << d << "v1b_v2a_std";

        outdata << std::endl;
    }

    for (int i = 0; i < p.numb_bin; i++)
    {

        double r = 0;

        if (p.lg_bins)
        {
            r = p.r_min * pow(10, (i + 0.5) * dlg_r);
        }
        else
        {
            r = p.r_min + (i + 0.5) * dr;
        }

        outdata << i;                        // 1
        outdata << d << r;                   // 2
        outdata << d << res.total.counts[i]; // 3

        if (p.r12_v1a)
            outdata << d << res.total.r12_v1a[i] << d << res.error.r12_v1a[i]; // 4,5
        if (p.r12_v1b)
            outdata << d << res.total.r12_v1b[i] << d << res.error.r12_v1b[i]; // 6,7
        if (p.r12_v2a)
            outdata << d << res.total.r12_v2a[i] << d << res.error.r12_v2a[i]; // 8,9
        if (p.r12_v2b)
            outdata << d << res.total.r12_v2b[i] << d << res.error.r12_v2b[i]; // 10,11

        if (p.v1a_v2a)
            outdata << d << res.total.v1a_v2a[i] << d << res.error.v1a_v2a[i]; // 12,13
        if (p.v1b_v2b)
            outdata << d << res.total.v1b_v2b[i] << d << res.error.v1b_v2b[i]; // 14,15
        if (p.v1a_v2b)
            outdata << d << res.total.v1a_v2b[i] << d << res.error.v1a_v2b[i]; // 16,17
        if (p.v1b_v2a)
            outdata << d << res.total.v1b_v2a[i] << d << res.error.v1b_v2a[i]; // 18,19

        if (i < p.numb_bin - 1)
        {
            outdata << std::endl;
        }
    }
}
