#include "SymFun.hpp"
#include <vector>
#include <cmath>

#ifdef DIM
#undef DIM
#endif
#define DIM 3

#ifdef MY_PI
#undef MY_PI
#endif
#define MY_PI 3.1415926535897932

int enzyme_dup;
int enzyme_out;
int enzyme_const;

inline double cut_cos(double const r, double const rcut) {
    return (r < rcut) ? 0.5 * (std::cos(MY_PI * r / rcut) + 1.0) : 0.0;
}
void sym_g2(double const eta, double const Rs, double const r, double const rcut, double &phi) {
    phi = std::exp(-eta * (r - Rs) * (r - Rs)) * cut_cos(r, rcut);
}

void symmetry_function_atomic(int const i,
                              double const *coords,
                              int const *particleSpeciesCodes,
                              int const *neighlist,
                              int const numnei,
                              double *const desc,
                              SymmetryFunctionParams *SymParam) {
    // prepare data
    VectorOfSizeDIM *coordinates = (VectorOfSizeDIM *) coords;
    int const iSpecies = particleSpeciesCodes[i];

    // Setup loop over neighbors of current particle
    for (int jj = 0; jj < numnei; ++jj) {
        // adjust index of particle neighbor
        int const j = neighlist[jj];
        int const jSpecies = particleSpeciesCodes[j];

        // cutoff between ij
        double rcutij = SymParam->rcut_2D_(iSpecies, jSpecies);

        // Compute rij
        double rij[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
            rij[dim] = coordinates[j][dim] - coordinates[i][dim];
        }

        double const rijsq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];
        double const rijmag = std::sqrt(rijsq);

        // if particles i and j not interact
        if (rijmag > rcutij) { continue; }

        // Loop over descriptors
        // two-body descriptors
        for (std::size_t p = 0; p < SymParam->name_.size(); ++p) {
            if (SymParam->name_[p] != 1 && SymParam->name_[p] != 2 && SymParam->name_[p] != 3) {
                continue;
            }

            int idx = SymParam->starting_index_[p];
            // Loop over same descriptor but different parameter set
            for (int q = 0; q < SymParam->num_param_sets_[p]; ++q) {
                double gc = 0.0;
                    double eta = SymParam->params_[p](q, 0);
                    auto Rs = SymParam->params_[p](q, 1);
                    sym_g2(eta, Rs, rijmag, rcutij, gc);
                desc[idx] += gc;
                ++idx;
            }
        }
    }
}

void grad_symmetry_function_atomic(int const i,
                                   double const *coords,
                                   double const *d_coords,
                                   int const *particleSpeciesCodes,
                                   int const *neighlist,
                                   int const numnei,
                                   double *const desc,
                                   double const *d_grad_loss_zeta,
                                   SymmetryFunctionParams *SymParam) {
    __enzyme_autodiff(symmetry_function_atomic,
                      enzyme_const, i,
                      enzyme_dup, coords, d_coords,
                      enzyme_const, particleSpeciesCodes,
                      enzyme_const, neighlist,
                      enzyme_const, numnei,
                      enzyme_dup, desc, d_grad_loss_zeta,
                      enzyme_const, SymParam);
}