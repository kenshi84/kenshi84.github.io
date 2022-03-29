#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;

#define MODE_STENCIL            // using a larger (2-ring neighborhood) stencil
//#define MODE_TWOSTEP            // two-step approach used in the GeoBrush paper

#define TEST_1D
//#define TEST_2D

const int N = 50;               // number of nodes
const int N_PLOT = 3;           // number of graph plots
const int N_ITER   = 10000;     // number of iterations for every plot (total number of iterations is N_PLOT * N_ITER)
const double DAMPING = 0.5;     // 0: no damping, 1: no update

struct Node {
    double value_;
    bool fixed_;
    Node()
        : value_(0)
        , fixed_(false)
    {}
    void fix(double value) {
        value_ = value;
        fixed_ = true;
    }
    operator double() const {
        return value_;
    }
};

#ifdef TEST_1D
double sample_func_1d(int i) {
    const double PI = 3 * acos(0.5);
    double x = i / (N - 1.);
    return sin(x * 2 * PI);
}

void test_1d() {
    vector<Node> f(N);
    f[0    ].fix(sample_func_1d(0    ));
    f[N - 1].fix(sample_func_1d(N - 1));
#if defined(MODE_STENCIL)
    f[1    ].fix(sample_func_1d(1    ));
    f[N - 2].fix(sample_func_1d(N - 2));
#elif defined(MODE_TWOSTEP)
    vector<Node> g(N - 1);
    g[0    ].fix(sample_func_1d(1    ) - sample_func_1d(0    ));
    g[N - 2].fix(sample_func_1d(N - 1) - sample_func_1d(N - 2));
#endif
    
    for (int p = 0; p < N_PLOT; ++p) {
        cout << "plot #: " << p << endl;
        for (int k = 0; k < N_ITER; ++k) {
            vector<Node> f_old = f;
#ifdef MODE_TWOSTEP
            vector<Node> g_old = g;
#endif
            for (int i = 0; i < N; ++i) {
                if (!f[i].fixed_) {
#if defined(MODE_STENCIL)
                    f[i].value_ = DAMPING * f_old[i] + (1 - DAMPING) / 6. * (
                        -     f_old[i - 2]
                        + 4 * f_old[i - 1]
                        + 4 * f_old[i + 1]
                        -     f_old[i + 2]);
#elif defined(MODE_TWOSTEP)
                    f[i].value_ = DAMPING * f_old[i] + (1 - DAMPING) * 0.5 * (
                        f_old[i - 1] + g_old[i - 1] +
                        f_old[i + 1] - g_old[i    ]);
#endif
                }
#ifdef MODE_TWOSTEP
                if (i < N - 1 && !g[i].fixed_) {
                    double g_p = g[i + 1].fixed_ ? g_old[i + 1] : (f_old[i + 2] - f_old[i + 1]);
                    double g_m = g[i - 1].fixed_ ? g_old[i - 1] : (f_old[i    ] - f_old[i - 1]);
                    g[i].value_ = DAMPING * g_old[i] + (1 - DAMPING) * 0.5 * (g_m + g_p);
                }
#endif
            }
        }
        stringstream ss;
#if defined(MODE_STENCIL)
        ss << "output_1d_stencil_" << p << ".txt";
#elif defined(MODE_TWOSTEP)
        ss << "output_1d_twostep_" << p << ".txt";
#endif
        ofstream ofs(ss.str().c_str());
        for (int i = 0; i < N; ++i) {
            ofs << f[i] << endl;
        }
    }
    /*
        to visualize the result, use the following command on gnuplot:
            plot "output_1d_stencil_0.txt" w l, "output_1d_stencil_1.txt" w l, "output_1d_stencil_2.txt" w l
            plot "output_1d_twostep_0.txt" w l, "output_1d_twostep_1.txt" w l, "output_1d_twostep_2.txt" w l
    */
}
#endif

#ifdef TEST_2D
struct Grid {
    int width_, height_;
    vector<Node> nodes_;
    Grid(int width, int height)
        : width_ (width)
        , height_(height)
        , nodes_ (width * height)
    {}
    const Node& operator()(int i, int j) const { return nodes_[width_ * j + i]; }
          Node& operator()(int i, int j)       { return nodes_[width_ * j + i]; }
};

double sample_func_2d(int i, int j) {
    double x = i / (N - 1.) - 0.5;
    double y = j / (N - 1.) - 0.5;
    return x * x + y * y;
}

void test_2d() {
    Grid f(N, N);
#ifdef MODE_TWOSTEP
    Grid gu(N - 1, N);
    Grid gv(N, N - 1);
#endif
    
    for (int k = 0; k < N; ++k) {
        f(k, 0    ).fix(sample_func_2d(k, 0    ));
        f(k, N - 1).fix(sample_func_2d(k, N - 1));
        f(0    , k).fix(sample_func_2d(0    , k));
        f(N - 1, k).fix(sample_func_2d(N - 1, k));
#ifdef MODE_STENCIL
        f(k, 1    ).fix(sample_func_2d(k, 1    ));
        f(k, N - 2).fix(sample_func_2d(k, N - 2));
        f(1    , k).fix(sample_func_2d(1    , k));
        f(N - 2, k).fix(sample_func_2d(N - 2, k));
#elif defined(MODE_TWOSTEP)
        gu(0    , k).fix(sample_func_2d(1    , k    ) - sample_func_2d(0    , k    ));
        gu(N - 2, k).fix(sample_func_2d(N - 1, k    ) - sample_func_2d(N - 2, k    ));
        gv(k, 0    ).fix(sample_func_2d(k    , 1    ) - sample_func_2d(k    , 0    ));
        gv(k, N - 2).fix(sample_func_2d(k    , N - 1) - sample_func_2d(k    , N - 2));
        if (k < N - 1) {
            gu(k, 0    ).fix(sample_func_2d(k + 1, 0    ) - sample_func_2d(k, 0    ));
            gu(k, N - 1).fix(sample_func_2d(k + 1, N - 1) - sample_func_2d(k, N - 1));
            gv(0    , k).fix(sample_func_2d(0    , k + 1) - sample_func_2d(0    , k));
            gv(N - 1, k).fix(sample_func_2d(N - 1, k + 1) - sample_func_2d(N - 1, k));
        }
#endif
    }
    
    for (int p = 0; p < N_PLOT; ++p) {
        cout << "plot #: " << p << endl;
        for (int k = 0; k < N_ITER; ++k) {
            Grid f_old = f;
#ifdef MODE_TWOSTEP
            Grid gu_old = gu;
            Grid gv_old = gv;
#endif
            for (int j = 0; j < N; ++j) for (int i = 0; i < N; ++i) {
                if (!f(i, j).fixed_) {
#if defined(MODE_STENCIL)
                    f(i, j).value_ = DAMPING * f_old(i, j) + (1 - DAMPING) / 20. * (
                         8 * (
                            f_old(i + 1, j) +
                            f_old(i - 1, j) +
                            f_old(i, j + 1) +
                            f_old(i, j - 1)) +
                        -2 * (
                            f_old(i + 1, j + 1) +
                            f_old(i - 1, j + 1) +
                            f_old(i + 1, j - 1) +
                            f_old(i - 1, j - 1)) +
                        -1 * (
                            f_old(i + 2, j) +
                            f_old(i - 2, j) +
                            f_old(i, j + 2) +
                            f_old(i, j - 2)));
#elif defined(MODE_TWOSTEP)
                    f(i, j).value_ = DAMPING * f_old(i, j) + (1 - DAMPING) * 0.25 * (
                        f_old(i - 1, j) + gu_old(i - 1, j) +
                        f_old(i + 1, j) - gu_old(i    , j) +
                        f_old(j, i - 1) + gv_old(j, i - 1) +
                        f_old(j, i + 1) - gv_old(j, i    ));
#endif
                }
#ifdef MODE_TWOSTEP
                if (i < N - 1 && !gu(i, j).fixed_) {
                    double gu_p0 = gu(i + 1, j).fixed_ ? gu_old(i + 1, j) : (f_old(i + 2, j    ) - f_old(i + 1, j    ));
                    double gu_m0 = gu(i - 1, j).fixed_ ? gu_old(i - 1, j) : (f_old(i    , j    ) - f_old(i - 1, j    ));
                    double gu_0p = gu(i, j + 1).fixed_ ? gu_old(i, j + 1) : (f_old(i + 1, j + 1) - f_old(i    , j + 1));
                    double gu_0m = gu(i, j - 1).fixed_ ? gu_old(i, j - 1) : (f_old(i + 1, j - 1) - f_old(i    , j - 1));
                    gu(i, j).value_ = DAMPING * gu_old(i, j) + (1 - DAMPING) * 0.25 * (gu_m0 + gu_p0 + gu_0p + gu_0m);
                }
                if (j < N - 1 && !gv(i, j).fixed_) {
                    double gv_p0 = gv(i + 1, j).fixed_ ? gv_old(i + 1, j) : (f_old(i + 1, j + 1) - f_old(i + 1, j    ));
                    double gv_m0 = gv(i - 1, j).fixed_ ? gv_old(i - 1, j) : (f_old(i - 1, j + 1) - f_old(i - 1, j    ));
                    double gv_0p = gv(i, j + 1).fixed_ ? gv_old(i, j + 1) : (f_old(i    , j + 2) - f_old(i    , j + 1));
                    double gv_0m = gv(i, j - 1).fixed_ ? gv_old(i, j - 1) : (f_old(i    , j    ) - f_old(i    , j - 1));
                    gv(i, j).value_ = DAMPING * gv_old(i, j) + (1 - DAMPING) * 0.25 * (gv_m0 + gv_p0 + gv_0p + gv_0m);
                }
#endif
            }
        }
        stringstream ss;
#if defined(MODE_STENCIL)
        ss << "output_2d_stencil_" << p << ".txt";
#elif defined(MODE_TWOSTEP)
        ss << "output_2d_twostep_" << p << ".txt";
#endif
        ofstream ofs(ss.str().c_str());
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i)
                ofs << f(i, j) << endl;
            ofs << endl;
        }
    }
    /*
        to visualize the result, use the following command on gnuplot:
            splot "output_2d_stencil_0.txt" w l, "output_2d_stencil_1.txt" w l, "output_2d_stencil_2.txt" w l
            splot "output_2d_twostep_0.txt" w l, "output_2d_twostep_1.txt" w l, "output_2d_twostep_2.txt" w l
    */
}
#endif

void main() {
#ifdef TEST_1D
    test_1d();
#endif
#ifdef TEST_2D
    test_2d();
#endif
}
