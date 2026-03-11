#ifndef PDMP_SUBSAMPLE_HPP
#define PDMP_SUBSAMPLE_HPP

#include <vector>
#include <stan/math/prim/fun/Eigen.hpp>

#ifdef _WIN32
#define PDMP_EXPORT __declspec(dllexport)
#else
#define PDMP_EXPORT __attribute__((visibility("default")))
#endif

namespace pdmp_subsample {
    static std::vector<int> indices_;
    static int m_ = 0;
}

extern "C" {
    PDMP_EXPORT void pdmp_set_subsample_indices(const int* idx, int new_m) {
        pdmp_subsample::indices_.assign(idx, idx + new_m);
        pdmp_subsample::m_ = new_m;
    }

    PDMP_EXPORT int pdmp_get_subsample_size() {
        return pdmp_subsample::m_;
    }
}

template <typename T>
inline Eigen::VectorXd get_subsampled_Y_real(const T& Y_full, std::ostream* pstream__) {
    int m = pdmp_subsample::m_;
    Eigen::VectorXd Y_sub(m);
    for (int i = 0; i < m; ++i) {
        Y_sub(i) = Y_full(pdmp_subsample::indices_[i]);
    }
    return Y_sub;
}

template <typename T>
inline std::vector<int> get_subsampled_Y_int(const T& Y_full, std::ostream* pstream__) {
    int m = pdmp_subsample::m_;
    std::vector<int> Y_sub(m);
    for (int i = 0; i < m; ++i) {
        Y_sub[i] = Y_full[pdmp_subsample::indices_[i]];
    }
    return Y_sub;
}

template <typename T>
inline Eigen::MatrixXd get_subsampled_Xc(const T& Xc_full, std::ostream* pstream__) {
    int m = pdmp_subsample::m_;
    int p = Xc_full.cols();
    Eigen::MatrixXd Xc_sub(m, p);
    for (int i = 0; i < m; ++i) {
        Xc_sub.row(i) = Xc_full.row(pdmp_subsample::indices_[i]);
    }
    return Xc_sub;
}

#endif
