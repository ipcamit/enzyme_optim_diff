#include <stdexcept>
#include <string>

#include <vector>
#include <cmath>
#define MY_PI 3.1415926535897932
#include <cstring>
#include <map>
#include <string>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <sstream>


template<class DataType = double>
class _Array_Basic {
public:
    _Array_Basic();
    _Array_Basic(std::size_t const count);
    _Array_Basic(std::size_t const count, DataType const value);
    _Array_Basic(std::size_t const count, DataType const *array);
    _Array_Basic(_Array_Basic<DataType> const &other);
    _Array_Basic(_Array_Basic<DataType> &&other);
    ~_Array_Basic();
    _Array_Basic<DataType> &operator=(_Array_Basic<DataType> const &other);
    _Array_Basic<DataType> &operator=(_Array_Basic<DataType> &&other);
    inline DataType const *data() const noexcept;
    inline DataType *data() noexcept;
    inline std::size_t size() const;
    inline void clear() noexcept;
    inline void shrink_to_fit();
    inline std::size_t capacity() const noexcept;
    inline void push_back(DataType const &value);
    inline void push_back(DataType &&value);
    std::vector<DataType> m;
protected:
    inline void _range_check(int _n) const;
    inline void _range_check(int _n, std::size_t tsize) const;
};

template<class DataType = double>
class Array1DView {
public:
    Array1DView(std::size_t const count, DataType *array);
    Array1DView(std::size_t const count, DataType const *array);
    Array1DView(Array1DView<DataType> const &other);
    ~Array1DView();
    inline DataType const *data() const noexcept;
    inline DataType *data() noexcept;
    inline const DataType operator()(int i) const;
    inline DataType &operator()(int i);
    inline DataType const at(int i) const;
    inline DataType &at(int i);
    const DataType operator[](int i) const;
    DataType &operator[](int i);
private:
    Array1DView() = delete;
    Array1DView<DataType> &operator=(Array1DView<DataType> const &other)
    = delete;
    Array1DView<DataType> &operator=(Array1DView<DataType> &&other) = delete;
protected:
    inline void _range_check(int _n, std::size_t tsize) const;
protected:
    /*! The extent of the container in the 1st mode */
    std::size_t _extentZero;
    /*! Data pointer */
    DataType *const m;
    };

template<class DataType = double>
class Array2DView {
public:
    Array2DView(std::size_t const extentZero, std::size_t const extentOne, DataType *array);
    Array2DView(std::size_t const extentZero, std::size_t const extentOne, DataType const *array);
    Array2DView(Array2DView<DataType> const &other);
    ~Array2DView();
    inline DataType const *data() const noexcept;
    inline DataType *data() noexcept;
    inline Array1DView<DataType> data_1D(int i);
    inline const DataType operator()(int i, int j) const;
    inline DataType &operator()(int i, int j);
    inline DataType const at(int i, int j) const;
    inline DataType &at(int i, int j);
    class j_operator {
    public:
        j_operator(Array2DView<DataType> &_array, int i);
        const DataType operator[](int j) const;
        DataType &operator[](int j);
    private:
        /*! Refernce to Array2D class */
        Array2DView<DataType> &j_array;
        std::size_t _i;
    };
    const j_operator operator[](int i) const;
    j_operator operator[](int i);
private:
    Array2DView() = delete;
    Array2DView<DataType> &operator=(Array2DView<DataType> const &other)
    = delete;
    Array2DView<DataType> &operator=(Array2DView<DataType> &&other) = delete;
protected:
    inline void _range_check(int _n, std::size_t tsize) const;
protected:
    /*! The extent of the container in the 1st mode */
    std::size_t _extentZero;
    /*! The extent of the container in the 2nd mode */
    std::size_t _extentOne;
    /*! Data pointer */
    DataType *const m;
};

template<class DataType = double>
class Array2D : public _Array_Basic<DataType> {
public:
    Array2D();
    Array2D(std::size_t const extentZero, std::size_t const extentOne);
    Array2D(std::size_t const extentZero, std::size_t const extentOne, DataType const value);
    Array2D(std::size_t const extentZero, std::size_t const extentOne, DataType const *array);
    Array2D(Array2D<DataType> const &other);
    Array2D(Array2D<DataType> &&other);
    ~Array2D();
    Array2D<DataType> &operator=(Array2D<DataType> const &other);
    Array2D<DataType> &operator=(Array2D<DataType> &&other);
    inline Array1DView<DataType> data_1D(int i);
    inline void resize(int const extentZero, int const extentOne);
    inline void resize(int const extentZero, int const extentOne, DataType const new_value);
    inline void resize(int const extentZero, int const extentOne, DataType const *new_array);
    inline const DataType operator()(int i, int j) const;
    inline DataType &operator()(int i, int j);
    inline DataType const at(int i, int j) const;
    inline DataType &at(int i, int j);
    class j_operator {
    public:
        j_operator(Array2D<DataType> &_array, int i);
        const DataType operator[](int j) const;
        DataType &operator[](int j);
    private:
        /*! Refernce to Array2D class */
        Array2D<DataType> &j_array;
        std::size_t _i;
};
    const j_operator operator[](int i) const;
    j_operator operator[](int i);
protected:
    /*! The extent of the container in the 1st mode */
    std::size_t _extentZero;
    /*! The extent of the container in the 2nd mode */
    std::size_t _extentOne;
};

void getNextDataLine(FILE *const filePtr, char *const nextLine, int const maxSize, int *endOfFileFlag);
int getXdouble(char *linePtr, const int N, double *list);
int getXint(char *linePtr, const int N, int *list);
void lowerCase(char *linePtr);
void lowerCase(std::string &InputLineArg);


#define DIM 3

typedef double VectorOfSizeDIM[DIM];


int enzyme_dup;
int enzyme_out;
int enzyme_const;


class SymmetryFunctionParams {
public:
    SymmetryFunctionParams();

    ~SymmetryFunctionParams();

    inline void set_species(std::vector<std::string> &species);

    inline void get_species(std::vector<std::string> &species);

    inline int get_num_species();
    int width;

    void set_cutoff(char const *name,
                    std::size_t Nspecies,
                    double const *rcut_2D);

    inline double get_cutoff(int iCode, int jCode);

    void add_descriptor(char const *name,
                        double const *values,
                        int row,
                        int col);

    int get_num_descriptors();

    std::vector<std::string> species_;
    std::vector<int> name_;
    std::vector<int> starting_index_;
    Array2D<double> rcut_2D_;
    std::vector<Array2D<double> > params_;
    std::vector<int> num_param_sets_;
    std::vector<int> num_params_;
    bool has_three_body_;
    double bhor2ang = 0.529177;
    void init();
};

inline double cut_cos(double r, double rcut);
void symmetry_function_atomic(int i,
                              double const *coords,
                              int const *particleSpeciesCodes,
                              int const *neighlist,
                              int numnei,
                              double * desc,
                              SymmetryFunctionParams *SymParam);

void sym_g2(double eta,
            double Rs,
            double r,
            double rcut,
            double &phi);


void __enzyme_autodiff(void (*)(
        int const,
        double const *,
        int const *,
        int const *,
        int const,
        double *const,
        SymmetryFunctionParams *),
                       int, int,
                       int, double const *, double const *,
                       int, int const *,
                       int, int const *,
                       int, int,
                       int, double *, double const *,
                       int, SymmetryFunctionParams *);

void grad_symmetry_function_atomic(int  i,
                                   double const *coords,
                                   double const *d_coords,
                                   int const *particleSpeciesCodes,
                                   int const *neighlist,
                                   int  numnei,
                                   double * desc,
                                   double const *d_grad_loss_zeta,
                                   SymmetryFunctionParams *SymParam);


SymmetryFunctionParams::SymmetryFunctionParams() : has_three_body_(false) {
}

SymmetryFunctionParams::~SymmetryFunctionParams() {}

inline void SymmetryFunctionParams::set_species(std::vector<std::string> &species) {
    species_.resize(species.size());
    std::copy(species.begin(), species.end(), species_.begin());
};

inline void SymmetryFunctionParams::get_species(std::vector<std::string> &species) {
    species.resize(species_.size());
    std::copy(species_.begin(), species_.end(), species.begin());
};

inline int SymmetryFunctionParams::get_num_species() { return species_.size(); }

void SymmetryFunctionParams::set_cutoff(char const *name,
                                        std::size_t const Nspecies,
                                        double const *rcut_2D) {
    (void) name;   // to avoid unused warning
    rcut_2D_.resize(Nspecies, Nspecies, rcut_2D);
}

inline double SymmetryFunctionParams::get_cutoff(int const iCode, int const jCode) {
    return rcut_2D_(iCode, jCode);
};

void SymmetryFunctionParams::add_descriptor(char const *name,
                                            double const *values,
                                            int const row,
                                            int const col) {
    // Enzyme string comparison workaround
    if (strcmp(name, "g1") == 0) { name_.push_back(1); };
    if (strcmp(name, "g2") == 0) { name_.push_back(2); };
    if (strcmp(name, "g3") == 0) { name_.push_back(3); };
    if (strcmp(name, "g4") == 0) { name_.push_back(4); };
    if (strcmp(name, "g5") == 0) { name_.push_back(5); };

    Array2D<double> params(row, col, values);
    params_.push_back(std::move(params));

    auto sum = std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
    starting_index_.push_back(sum);

    num_param_sets_.push_back(row);
    num_params_.push_back(col);

    if (strcmp(name, "g4") == 0 || strcmp(name, "g5") == 0) {
        has_three_body_ = true;
    }
}

int SymmetryFunctionParams::get_num_descriptors() {
    return std::accumulate(num_param_sets_.begin(), num_param_sets_.end(), 0);
}


void SymmetryFunctionParams::init() {
    int n_species = 1;
    double * cutoff_matrix = new double[n_species * n_species];
    for (int i=0;i<n_species;i++){
        for (int j=0;j<n_species;j++){
            *( cutoff_matrix + n_species * i + j) = 3.77;
        }
    }
    std::string cutoff_function = "cos";
    width = 51;
    int n_func = 2;

    std::vector<std::string> sym_func_list = {"g2", "g4"};

    std::vector<int> sym_func_lengths = {16, 129};

    std::vector<std::vector<double>> sym_func_elements = {{0.001, 0.0, 0.01, 0.0, 0.02, 0.0, 0.035, 0.0, 0.06, 0.0, 0.1, 0.0, 0.2, 0.0, 0.4, 0.0},
    {1., -1., 0.0001, 1., 1., 0.0001, 2., -1., 0.0001, 2., 1., 0.0001, 1., -1., 0.003, 1., 1., 0.003, 2, -1, 0.003, 2, 1, 0.003, 1, -1, 0.008, 1, 1, 0.008, 2, -1, 0.008, 2, 1, 0.008, 1, -1, 0.015, 1, 1, 0.015, 2, -1, 0.015, 2, 1, 0.015, 4, -1, 0.015, 4, 1, 0.015, 16, -1, 0.015, 16, 1, 0.015, 1, -1, 0.025, 1, 1, 0.025, 2, -1, 0.025, 2, 1, 0.025, 4, -1, 0.025, 4, 1, 0.025, 16, -1, 0.025, 16, 1, 0.025, 1, -1, 0.045, 1, 1, 0.045, 2, -1, 0.045, 2, 1, 0.045, 4, -1, 0.045, 4, 1, 0.045, 16, -1, 0.045, 16, 1, 0.045, 1, -1, 0.08, 1, 1, 0.08, 2, -1, 0.08, 2, 1, 0.08, 4, -1, 0.08, 4, 1, 0.08, 16, 1, 0.08} };



    for (int i = 0; i < n_func ; i++){
        if (sym_func_list[i] == "g2"){
            for (int j = 0; j < sym_func_elements[i].size(); j = j+2) sym_func_elements[i][j] /= (bhor2ang * bhor2ang);
        } else if (sym_func_list[i] == "g4") {
            for (int j = 2; j < sym_func_elements[i].size(); j = j+3) sym_func_elements[i][j] /= (bhor2ang * bhor2ang);
        }
    }

    std::vector<int> dims = {8, 2, 43, 3};
    set_cutoff(cutoff_function.c_str(), n_species, cutoff_matrix);

    for (int i =0 ; i< n_func; i++){
        add_descriptor(sym_func_list[i].c_str(), sym_func_elements[i].data(), dims[2 *i], dims[2*i+1] );
    }
    delete[] cutoff_matrix;
}


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

int main() {

    int neigh_list[] = {120,  70,  15,  80,  155,  2,  153, 136,  144,  6,  114,  111,  137,  138,  5};
    int n_neigh[] = {7, 8};
    double desc[] = {2.38468, 1.99139, 1.63009, 1.20736, 0.732224, 0.329115, 0.0446837, 0.000832121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.49234, 2.09188, 1.72222, 1.28688, 0.792281, 0.365058, 0.0528919, 0.00113331, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double d_desc[] = {-0.0459296, -0.0571872, -0.0329145, -0.0385561, -0.0416095, -0.0497505, -0.0460529, -0.0499233, -0.0223421, -0.0228727, -0.0252878, -0.0316037, -0.0222051, -0.0243794, -0.0250438, -0.0121398, -0.0215549, -0.0122885, -0.0155702, -0.0216616, -0.00928747, -0.0239285, -0.0203377, -0.0201377, -0.0272253, -0.026079, -0.014359, -0.022687, -0.0222708, -0.0300837, -0.0266789, -0.0303648, -0.0289803, -0.0291265, -0.0127385, -0.0204652, -0.0258385, -0.0396512, -0.0178335, -0.0225855, -0.0262925, -0.0188549, 0.00172229, -0.0242082, -0.0194935, -0.0189299, -0.0226348, -0.0205645, -0.0181264, -0.0277141, -0.0157652, -0.0388413, -0.0482813, -0.0277865, -0.032564, -0.0350948, -0.0420019, -0.0388929, -0.0421405, -0.0189227, -0.0193064, -0.021368, -0.0267672, -0.018772, -0.0206341, -0.0211213, -0.0102563, -0.0182032, -0.0103453, -0.0131602, -0.018283, -0.00783743, -0.0201724, -0.017121, -0.0170211, -0.0229717, -0.0219972, -0.0121846, -0.0191228, -0.0187706, -0.0254743, -0.0225702, -0.0256507, -0.0244421, -0.0246517, -0.0108475, -0.0172759, -0.0217582, -0.0334613, -0.0150477, -0.0190535, -0.0222183, -0.0159371, 0.00147253, -0.0204223, -0.0164294, -0.0160068, -0.0191025, -0.0173686, -0.0152722, -0.0233994, -0.0133275};

    double coords[] = {0.244549, 0.1286, 10.6624, 1.73029, 1.36074, 1.17815, 2.95995, 2.72051, 10.699, 6.90013, 4.00404, 3.95301, 2.82161, 10.7364, 2.55305, 3.93608, 1.13148, 4.10265, 0.279865, 2.79618, 2.43383, 1.72218, 4.08926, 3.86036, 0.209679, 0.0136685, 5.27347, 1.48183, 1.34708, 6.66339, 2.81651, 2.7307, 5.36199, 4.23699, 4.04211, 6.63741, 2.96636, 10.8597, 8.03465, 4.29722, 1.42171, 9.32751, 9.73274, 4.07308, 6.51836, 0.0370598, 2.65909, 7.94201, 10.8549, 5.36131, 10.7369, 1.39826, 6.7256, 1.24454, 2.78667, 8.1682, 10.7896, 4.2552, 9.45391, 1.2355, 2.95533, 5.46094, 2.53352, 4.31299, 6.82168, 3.88275, 0.0623023, 8.17698, 2.58624, 1.50293, 9.41922, 3.93406, 0.233928, 5.44444, 5.1516, 1.57376, 6.73965, 6.62621, 2.98098, 8.15044, 5.34163, 3.5635, 10.3251, 5.85504, 2.77497, 5.28587, 7.96171, 3.94244, 6.6129, 9.49359, 0.229244, 8.12558, 7.99351, 1.59443, 9.50522, 9.3491, 5.53332, 10.8364, 10.6917, 5.56799, 2.80038, 2.47042, 6.90624, 1.33905, 1.26044, 9.56892, 3.95367, 1.27585, 9.76477, 1.51457, 3.78284, 7.02589, 1.38808, 6.65603, 1.39344, 3.97548, 9.38958, 4.27015, 4.20669, 1.06284, 6.21818, 0.446436, 4.74922, 9.61209, 1.33158, 9.36732, 8.3117, 2.71327, 5.20123, 7.05119, 4.06019, 9.2623, 8.29674, 10.7298, 8.05205, 9.69008, 9.67516, 1.20536, 5.63507, 2.76046, 7.99448, 8.23299, 2.63751, 10.7358, 5.77263, 5.46, 10.6353, 7.04579, 6.82197, 1.19108, 8.37253, 8.2257, 10.6973, 8.33411, 0.123501, 2.55475, 8.21039, 5.34252, 2.62742, 9.65035, 6.72057, 3.87815, 5.65865, 8.09068, 2.5227, 7.03813, 9.43229, 3.87684, 5.56006, 5.45311, 5.29477, 6.95069, 6.77211, 6.69254, 8.32697, 8.07264, 5.3524, 9.67322, 9.46285, 6.65913, 8.40108, 5.44149, 7.93652, 9.74047, 6.80838, 9.35191, 5.72574, 7.66117, 8.47106, 6.92554, 9.43493, 9.35384, -2.48947, -2.6363, -0.164651, -2.56526, -0.132152, 8.05205, -1.17192, -1.18684, 1.20536, -2.48947, -2.6363, 10.6973, -2.53503, -2.78936, 5.3524, -1.18878, -1.39915, 6.65913, -1.17192, -1.18684, 12.0674, -0.00706323, 5.36131, -0.12505, -1.24991, 1.33158, -1.49468, -2.62901, 2.63751, -0.126211, -2.48947, 8.2257, -0.164651, -1.12153, 6.80838, -1.51009, -1.12926, 4.07308, 6.51836, -0.00706323, 5.36131, 10.7369, -1.29308, 3.95367, 1.27585, -1.09723, 1.51457, 3.78284, -1.24991, 1.33158, 9.36732, -2.5503, 2.71327, 5.20123, -2.56526, 10.7298, 8.05205, -1.17192, 9.67516, 1.20536, -2.62901, 2.63751, 10.7358, -2.48947, 8.2257, 10.6973, -2.52789, 0.123501, 2.55475, -2.65161, 5.34252, 2.62742, -1.21165, 6.72057, 3.87815, -2.53503, 8.07264, 5.3524, -1.18878, 9.46285, 6.65913, -2.46092, 5.44149, 7.93652, -1.12153, 6.80838, 9.35191, -1.29308, 3.95367, 12.1378, -1.17192, 9.67516, 12.0674, -2.52789, 0.123501, 13.4168, -2.65161, 5.34252, 13.4894, -1.24991, 12.1936, -1.49468, -2.62901, 13.4995, -0.126211, -1.09723, 12.3766, 3.78284, -1.24991, 12.1936, 9.36732, -2.5503, 13.5753, 5.20123, -2.62901, 13.4995, 10.7358, -2.52789, 10.9855, 2.55475, -2.52789, 10.9855, 13.4168, 2.78667, -2.6938, -0.072409, 1.59443, -1.35678, -1.5129, 5.53332, -0.0256437, -0.170315, 8.37253, -2.6363, -0.164651, 5.72574, -3.20083, -2.39094, 6.92554, -1.42707, -1.50816, 2.82161, -0.125615, 2.55305, 2.96636, -0.00230579, 8.03465, 2.78667, -2.6938, 10.7896, 4.2552, -1.40809, 1.2355, 0.0623023, -2.68502, 2.58624, 1.50293, -1.44278, 3.93406, 2.98098, -2.71156, 5.34163, 3.5635, -0.53693, 5.85504, 0.229244, -2.73642, 7.99351, 1.59443, -1.35678, 9.3491, 5.53332, -0.0256437, 10.6917, 8.29674, -0.132152, 8.05205, 9.69008, -1.18684, 1.20536, 8.37253, -2.6363, 10.6973, 5.65865, -2.77132, 2.5227, 7.03813, -1.42971, 3.87684, 8.32697, -2.78936, 5.3524, 9.67322, -1.39915, 6.65913, 5.72574, -3.20083, 8.47106, 6.92554, -1.42707, 9.35384, 2.82161, -0.125615, 13.415, 4.2552, -1.40809, 12.0975, 0.0623023, -2.68502, 13.4482, 9.69008, -1.18684, 12.0674, 5.65865, -2.77132, 13.3847, 0.244549, 0.1286, -0.199593, 2.95995, 2.72051, -0.162971, 4.29722, 1.42171, -1.53449, 10.8549, 5.36131, -0.12505, 2.78667, 8.1682, -0.072409, 3.94244, 6.6129, -1.36841, 1.59443, 9.50522, -1.5129, 5.53332, 10.8364, -0.170315, 1.39344, 3.97548, -1.47242, 9.61209, 1.33158, -1.49468, 7.05119, 4.06019, -1.5997, 8.23299, 2.63751, -0.126211, 5.77263, 5.46, -0.22672, 8.37253, 8.2257, -0.164651, 9.74047, 6.80838, -1.51009, 5.72574, 7.66117, -2.39094, 6.92554, 9.43493, -1.50816, 1.73029, 1.36074, 12.0401, 2.82161, 10.7364, 13.415, 0.279865, 2.79618, 13.2958, 1.39826, 6.7256, 12.1065, 4.2552, 9.45391, 12.0975, 2.95533, 5.46094, 13.3955, 0.0623023, 8.17698, 13.4482, 5.56799, 2.80038, 13.3324, 6.90624, 1.33905, 12.1224, 9.56892, 3.95367, 12.1378, 4.27015, 4.20669, 11.9248, 9.69008, 9.67516, 12.0674, 7.04579, 6.82197, 12.0531, 8.33411, 0.123501, 13.4168, 8.21039, 5.34252, 13.4894, 5.65865, 8.09068, 13.3847, 0.244549, 10.9906, -0.199593, 2.95995, 13.5825, -0.162971, 4.29722, 12.2837, -1.53449, 9.61209, 12.1936, -1.49468, 8.23299, 13.4995, -0.126211, 0.244549, 10.9906, 10.6624, 1.73029, 12.2227, 1.17815, 2.95995, 13.5825, 10.699, 3.93608, 11.9935, 4.10265, 0.279865, 13.6582, 2.43383, 0.209679, 10.8757, 5.27347, 1.48183, 12.2091, 6.66339, 2.81651, 13.5927, 5.36199, 4.29722, 12.2837, 9.32751, 0.0370598, 13.5211, 7.94201, 5.56799, 13.6624, 2.47042, 6.90624, 12.201, 1.26044, 9.76477, 12.3766, 3.78284, 7.02589, 12.2501, 6.65603, 6.21818, 11.3084, 4.74922, 9.61209, 12.1936, 9.36732, 8.3117, 13.5753, 5.20123, 5.63507, 13.6225, 7.99448, 8.23299, 13.4995, 10.7358, 8.33411, 10.9855, 2.55475, 1.73029, 12.2227, 12.0401, 0.279865, 13.6582, 13.2958, 5.56799, 13.6624, 13.3324, 6.90624, 12.201, 12.1224, 8.33411, 10.9855, 13.4168, 13.6487, -2.6938, -0.072409, 12.4564, -1.35678, -1.5129, 13.6836, -0.125615, 2.55305, 13.8284, -0.00230579, 8.03465, 13.6487, -2.6938, 10.7896, 10.9243, -2.68502, 2.58624, 12.3649, -1.44278, 3.93406, 13.843, -2.71156, 5.34163, 14.4255, -0.53693, 5.85504, 11.0912, -2.73642, 7.99351, 12.4564, -1.35678, 9.3491, 13.6836, -0.125615, 13.415, 10.9243, -2.68502, 13.4482, 11.1065, 0.1286, -0.199593, 13.822, 2.72051, -0.162971, 13.6487, 8.1682, -0.072409, 12.4564, 9.50522, -1.5129, 12.2554, 3.97548, -1.47242, 11.1065, 0.1286, 10.6624, 12.5923, 1.36074, 1.17815, 13.822, 2.72051, 10.699, 13.6836, 10.7364, 2.55305, 11.1419, 2.79618, 2.43383, 12.5842, 4.08926, 3.86036, 11.0717, 0.0136685, 5.27347, 12.3438, 1.34708, 6.66339, 13.6785, 2.7307, 5.36199, 13.8284, 10.8597, 8.03465, 10.8991, 2.65909, 7.94201, 12.2603, 6.7256, 1.24454, 13.6487, 8.1682, 10.7896, 13.8173, 5.46094, 2.53352, 10.9243, 8.17698, 2.58624, 12.3649, 9.41922, 3.93406, 11.0959, 5.44444, 5.1516, 12.4358, 6.73965, 6.62621, 13.843, 8.15044, 5.34163, 14.4255, 10.3251, 5.85504, 13.637, 5.28587, 7.96171, 11.0912, 8.12558, 7.99351, 12.4564, 9.50522, 9.3491, 12.2554, 3.97548, 9.38958, 12.5923, 1.36074, 12.0401, 13.6836, 10.7364, 13.415, 11.1419, 2.79618, 13.2958, 12.2603, 6.7256, 12.1065, 13.8173, 5.46094, 13.3955, 10.9243, 8.17698, 13.4482, 11.1065, 10.9906, -0.199593, 13.822, 13.5825, -0.162971, 11.1065, 10.9906, 10.6624, 12.5923, 12.2227, 1.17815, 13.822, 13.5825, 10.699, 11.1419, 13.6582, 2.43383, 11.0717, 10.8757, 5.27347, 12.3438, 12.2091, 6.66339, 13.6785, 13.5927, 5.36199, 10.8991, 13.5211, 7.94201, 12.5923, 12.2227, 12.0401, 11.1419, 13.6582, 13.2958};

    double forces[259 * 3];
    int species[259];
    int all_data = 777;
    for(int i =0;i<777;i++) forces[i] = 0.;
    for(int i =0;i<259;i++) species[i] = 0;
    int neigh_from = 0;

    auto sf = SymmetryFunctionParams();
    sf.init();

    for(int i = 0; i < 2; i++) {
        auto n_neigh_ = n_neigh[i];
        std::vector<int> n_list(neigh_list + neigh_from,
                                            neigh_list + neigh_from + n_neigh_);
        neigh_from += n_neigh_;
        grad_symmetry_function_atomic(i,
                                    coords,
                                    forces,
                                    species,
                                    n_list.data(),
                                    n_neigh_,
                                    desc + (i * 51),
                                    d_desc + (i * 51),
                                    &sf);
    }

    for (int i = 0; i < 3; i++){
        std::cout << forces[3 * i + 0] << "  " << forces[3 * i + 1] << "  " << forces[3 * i + 2] <<"\n";
    }
    return 0;
}



// --------------------------- Implementation --------------------------- //

template<class DataType>
_Array_Basic<DataType>::_Array_Basic() {}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count) :
        m(count, static_cast<DataType>(0)) {}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count,
                                     DataType const value) :
        m(count, value) {}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count,
                                     DataType const *array) :
        m(array, array + count) {}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> const &other) :
        m(other.m) {}

template<class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> &&other) :
        m(std::move(other.m)) {}

template<class DataType>
_Array_Basic<DataType>::~_Array_Basic() {}

template<class DataType>
_Array_Basic<DataType> &
_Array_Basic<DataType>::operator=(_Array_Basic<DataType> const &other) {
    m.resize(other.size());
    std::copy(other.m.begin(), other.m.end(), m.begin());
    return *this;}

template<class DataType>
_Array_Basic<DataType> &
_Array_Basic<DataType>::operator=(_Array_Basic<DataType> &&other) {
    m = std::move(other.m);
    return *this;}

template<class DataType>
inline DataType const *_Array_Basic<DataType>::data() const noexcept {
    return m.data();}

template<class DataType>
inline DataType *_Array_Basic<DataType>::data() noexcept {
    return m.data();}

template<class DataType>
inline std::size_t _Array_Basic<DataType>::size() const {
    return m.size();}

template<class DataType>
inline void _Array_Basic<DataType>::clear() noexcept {
    m.clear();}

template<class DataType>
inline void _Array_Basic<DataType>::shrink_to_fit() {
    m.shrink_to_fit();}

template<class DataType>
inline std::size_t _Array_Basic<DataType>::capacity() const noexcept {
    return m.capacity();}

template<class DataType>
inline void _Array_Basic<DataType>::push_back(DataType const &value) {
    m.push_back(value);}

template<class DataType>
inline void _Array_Basic<DataType>::push_back(DataType &&value) {
    m.push_back(value);}

template<class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n) const {
    if (_n >= size()) {
        std::abort();
    }}

template<class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n,
                                                 std::size_t tsize) const {
    if (_n >= tsize) {
        std::abort();
    }}

template<class DataType>
Array1DView<DataType>::Array1DView(std::size_t const count, DataType *array) :
        _extentZero(count), m(array) {}

template<class DataType>
Array1DView<DataType>::Array1DView(Array1DView<DataType> const &other) :
        _extentZero(other._extentZero), m(other.m) {}

template<class DataType>
Array1DView<DataType>::~Array1DView() {}

template<class DataType>
inline DataType const *Array1DView<DataType>::data() const noexcept {
    return m;}

template<class DataType>
inline DataType *Array1DView<DataType>::data() noexcept {
    return m;}

template<class DataType>
inline const DataType Array1DView<DataType>::operator()(int i) const {
    return m[i];}

template<class DataType>
inline DataType &Array1DView<DataType>::operator()(int i) {
    return m[i];}

template<class DataType>
inline DataType &Array1DView<DataType>::at(int i) {
    _range_check(i, _extentZero);
    return m[i];}

template<class DataType>
inline DataType const Array1DView<DataType>::at(int i) const {
    _range_check(i, _extentZero);
    return m[i];}

template<class DataType>
const DataType Array1DView<DataType>::operator[](int i) const {
    return m[i];}

template<class DataType>
DataType &Array1DView<DataType>::operator[](int i) {
    return m[i];}

template<class DataType>
inline void Array1DView<DataType>::_range_check(int _n, std::size_t tsize) const {
    if (_n >= tsize) {
        std::abort();
    }
}

template<class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType *array) :
        _extentZero(extentZero), _extentOne(extentOne), m(array) {}

template<class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType const *array) :
        _extentZero(extentZero),
        _extentOne(extentOne),
        m(const_cast<DataType *>(array)) {}

template<class DataType>
Array2DView<DataType>::Array2DView(Array2DView<DataType> const &other) :
        _extentZero(other._extentZero), _extentOne(other._extentOne), m(other.m) {}

template<class DataType>
Array2DView<DataType>::~Array2DView() {}

template<class DataType>
inline DataType const *Array2DView<DataType>::data() const noexcept {
    return m;}

template<class DataType>
inline DataType *Array2DView<DataType>::data() noexcept {
    return m;}

template<class DataType>
inline Array1DView<DataType> Array2DView<DataType>::data_1D(int i) {
    return Array1DView<DataType>(_extentOne, m + i * _extentOne);}

template<class DataType>
inline const DataType Array2DView<DataType>::operator()(int i, int j) const {
    std::size_t const _n = i * _extentOne + j;
    return m[_n];}

template<class DataType>
inline DataType &Array2DView<DataType>::operator()(int i, int j) {
    std::size_t const _n = i * _extentOne + j;
    return m[_n];}

template<class DataType>
inline DataType &Array2DView<DataType>::at(int i, int j) {
    _range_check(i, _extentZero);
    _range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return m[_n];}

template<class DataType>
inline DataType const Array2DView<DataType>::at(int i, int j) const {
    _range_check(i, _extentZero);
    _range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return m[_n];}

template<class DataType>
Array2DView<DataType>::j_operator::j_operator(Array2DView<DataType> &_array,
                                              int i) :
        j_array(_array), _i(i) {}

template<class DataType>
const DataType Array2DView<DataType>::j_operator::operator[](int j) const {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];}

template<class DataType>
DataType &Array2DView<DataType>::j_operator::operator[](int j) {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];}

template<class DataType>
const typename Array2DView<DataType>::j_operator
Array2DView<DataType>::operator[](int i) const {
    return j_operator(*this, i);}

template<class DataType>
typename Array2DView<DataType>::j_operator
Array2DView<DataType>::operator[](int i) {
    return j_operator(*this, i);}

template<class DataType>
inline void Array2DView<DataType>::_range_check(int _n, std::size_t tsize) const {
    if (_n >= tsize) {
        std::abort();
}}

template<class DataType>
Array2D<DataType>::Array2D() :
        _Array_Basic<DataType>(), _extentZero(0), _extentOne(0) {}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne) :
        _Array_Basic<DataType>(extentZero * extentOne),
        _extentZero(extentZero),
        _extentOne(extentOne) {}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const value) :
        _Array_Basic<DataType>(extentZero * extentOne, value),
        _extentZero(extentZero),
        _extentOne(extentOne) {}

template<class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const *array) :
        _Array_Basic<DataType>(extentZero * extentOne, array),
        _extentZero(extentZero),
        _extentOne(extentOne) {}

template<class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> const &other) :
        _Array_Basic<DataType>(other),
        _extentZero(other._extentZero),
        _extentOne(other._extentOne) {}

template<class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> &&other) :
        _Array_Basic<DataType>(std::move(other)),
        _extentZero(other._extentZero),
        _extentOne(other._extentOne) {}

template<class DataType>
Array2D<DataType>::~Array2D() {}

template<class DataType>
Array2D<DataType> &
Array2D<DataType>::operator=(Array2D<DataType> const &other) {
    _Array_Basic<DataType>::operator=(other);
    _extentZero = other._extentZero;
    _extentOne = other._extentOne;
    return *this;}

template<class DataType>
Array2D<DataType> &Array2D<DataType>::operator=(Array2D<DataType> &&other) {
    _Array_Basic<DataType>::operator=(std::move(other));
    _extentZero = other._extentZero;
    _extentOne = other._extentOne;
    return *this;}

template<class DataType>
inline Array1DView<DataType> Array2D<DataType>::data_1D(int i) {
    return Array1DView<DataType>(_extentOne, this->m.data() + i * _extentOne);}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero, int const extentOne) {
    _extentZero = extentZero;
    _extentOne = extentOne;
    std::size_t const _n = _extentZero * _extentOne;
    this->m.resize(_n, static_cast<DataType>(0));
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const new_value) {
    _extentZero = extentZero;
    _extentOne = extentOne;
    std::size_t const _n = _extentZero * _extentOne;
    this->m.resize(_n, new_value);
}

template<class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const *new_array) {
    _extentZero = extentZero;
    _extentOne = extentOne;
    std::size_t const _n = _extentZero * _extentOne;
    this->m.resize(_n);
    std::copy(new_array, new_array + _n, this->m.data());
}

template<class DataType>
inline const DataType Array2D<DataType>::operator()(int i, int j) const {
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
inline DataType &Array2D<DataType>::operator()(int i, int j) {
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
inline DataType &Array2D<DataType>::at(int i, int j) {
    this->_range_check(i, _extentZero);
    this->_range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
inline DataType const Array2D<DataType>::at(int i, int j) const {
    this->_range_check(i, _extentZero);
    this->_range_check(j, _extentOne);
    std::size_t const _n = i * _extentOne + j;
    return this->m[_n];
}

template<class DataType>
Array2D<DataType>::j_operator::j_operator(Array2D<DataType> &_array, int i) :
        j_array(_array), _i(i) {
}

template<class DataType>
const DataType Array2D<DataType>::j_operator::operator[](int j) const {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];
}

template<class DataType>
DataType &Array2D<DataType>::j_operator::operator[](int j) {
    std::size_t const _n = _i * j_array._extentOne + j;
    return j_array.m[_n];
}

template<class DataType>
const typename Array2D<DataType>::j_operator
Array2D<DataType>::operator[](int i) const {
    return j_operator(*this, i);
}

template<class DataType>
typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i) {
    return j_operator(*this, i);
}

