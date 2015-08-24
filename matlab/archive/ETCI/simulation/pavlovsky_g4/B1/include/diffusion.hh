#ifndef DIFFUSION_HH
#define DIFFUSION_HH
#include "math.h"
#include "TH2D.h"
#include "TRandom.h"

class diffusion
{
public:
    diffusion(TH2D* hist, double pix_size_um);
    void diffuseInteraction(double x0_um, double y0_um, double z0_um, double Edep_eV);

    double getDriftTime(){ return drift_time; }
    double getSatParam(){ return sat_param; }
    double pointSpreadFunction(double x, double y, void *data);
    double gauss_legendre_2D_cube(int n, void *data, double a, double b, double c, double d);

private:
    void driftTime(double zeta, int alpha = 1);
    void checkCoordinates(double x, double y, double z);
    double saturationParameter(double bias_V = 200., double gamma = 1.5 );

    // function that returns point density at (x, y)

    TRandom *prng;
    double eps_Si_F_cm;
    double N_a_cm3;
    double drift_time;
    double ri2;
    double temp_k;
    double mu0_um2;
    double diff_const;
    double ccd_depth_um;
    double pos_um[3];
    double q_0;
    double to;
    double sat_param;
    double pixel_size_um;
    double rd2;
    TH2D* ccdHist;


};

#endif // DIFFUSION_HH
