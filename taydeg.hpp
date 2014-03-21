#include <algorithm>

/**
 * Computes the degree N for the Taylor polynomial
 * of exp(tP) to have error less than c*eps*exp(t),
 * where c<1/2, then returns eps
 *
 * ( so exp(-t(I-P)) has error less than eps )
 */
unsigned int taylordegree(const double t, double& eps) {
    double eps_exp_t = eps*exp(t);
    double error = exp(t)-1.;
    double last = 1.;
    double k = 0.;
    while(error > eps_exp_t/2){
        k = k + 1.;
        last = (last*t)/k;
        error = error - last;
    }
    eps = (eps_exp_t - error)/exp(t);
    return std::max((int)k, (int)1);
}
