#include "diffusion.hh"
#include "QDebug"
#include "tr1/random"

diffusion::diffusion(TH2D* hist, double pix_size_um)
{
    prng = new TRandom(3243);
    ccdHist = hist;
    eps_Si_F_cm = 8.85E-14*11.68;
    N_a_cm3 = 1.3E13;
    drift_time = -1;
    ri2 = 1.;
    temp_k = 130.;
    mu0_um2 = 2714.2634322713E8; // from mobility.ods/xls
    diff_const = mu0_um2*(0.0258*temp_k/300.);
    ccd_depth_um = 650.0;
    pos_um[3] = 0;
    q_0 = 0;
    to = eps_Si_F_cm/(1.602E-19*N_a_cm3*mu0_um2/1.E8);
    sat_param = saturationParameter();
    pixel_size_um = pix_size_um;

}

void diffusion::diffuseInteraction(double x0_um, double y0_um, double z0_um, double Edep_eV)
{
    // eV/e-h Temperature dependence of the average energy per electron-hole pair in Silicon and Germanium
    // R. Stuck , J. P. Ponpon , R. Berger & P. Siffert
    q_0 = /*prng->PoissonD(*/Edep_eV/3.705/*)*/;

    ri2 = 0.000144 * pow(q_0*3.705/1000., 3.5);

    //    for(double x = pos_um[0]-10.*pixel_size_um; x < pos_um[0]+10.*pixel_size_um; x+=0.5){
    //        for( double y=pos_um[1]-10.*pixel_size_um; y<pos_um[1]+10.*pixel_size_um; y+=0.5){
    // Transform to correct coordinates
    checkCoordinates(z0_um,y0_um,x0_um);
    driftTime(x0_um/ccd_depth_um,1);
    if(0.==drift_time) {
        qDebug() << "invalid zeta no diffusion";
        return;
    }
    double w=0, sum=0;
    double window_width = 4.*sqrt(rd2+ri2)/sqrt(2.);
    double grid = 0.125;
    for(double y = pos_um[0]-window_width; y < pos_um[1]+window_width; y+=grid){
        for( double x=pos_um[1]-window_width; x<pos_um[0]+window_width; x+=grid){
            w=( pointSpreadFunction(x,y,NULL)*grid*grid );
            if( !isnan(w) && !isinf(w)){ //figure out why these are misbehaving
                sum+=w;
                ccdHist->Fill(x/pixel_size_um, y/pixel_size_um,w*3.705);
            }
        }
    }
//    qDebug() << 2*sqrt(rd2+ri2)/sqrt(2.);
//    qDebug() << Edep_eV << diff_const << drift_time << rd2 << ri2;
}

// Charge diffusion in CCD X-ray detectors
// George G. Pavlov, John A. Nousek; Penn State
void diffusion::driftTime(double zeta, int alpha)
{
    if(1 == alpha){
        // if 1-zeta is >= 0 then the hit is outside depletion, throw away
        if( 0 >= 1-zeta ){ return; }
        drift_time = to*( sat_param*zeta - log(1-zeta) );
        rd2 = 4*diff_const*drift_time;
    }else{
        //implement alpha for other values
    }
}

// Coordinate system is different for diffusion with the x, z variables
// permuted. Sorry about this but fixed in function call in diffuseInteraction
void diffusion::checkCoordinates(double x, double y, double z)
{
    pos_um[0] = x+726.*pixel_size_um/2;
    pos_um[1] = y+1454.*pixel_size_um/2;
    pos_um[2] = z;
}

double diffusion::saturationParameter(double bias_V, double gamma)
{
//    return 2.5*sqrt(N_a_cm3/1.E13)*sqrt(bias_V/10.)/pow(temp_k/153, gamma);
    return 1.602E-19*N_a_cm3*ccd_depth_um/(1.E4)/eps_Si_F_cm/(2.E3*pow(temp_k/153.,1.5));
}

double diffusion::pointSpreadFunction(double x, double y, void* data)
{
    double r2 = pow(x-pos_um[0], 2) + pow(y-pos_um[1],2);
    return  q_0/(M_PI*(ri2+rd2))*exp( - ( r2 ) / (ri2+rd2) );
}

double diffusion::gauss_legendre_2D_cube(int n, void *data, double a, double b, double c, double d)
{
////    /* n = 16 */
//    static double x[8] = {0.0950125098376374401853193,0.2816035507792589132304605,0.4580167776572273863424194,0.6178762444026437484466718,0.7554044083550030338951012,0.8656312023878317438804679,0.9445750230732325760779884,0.9894009349916499325961542};
//    static double w[8] = {0.1894506104550684962853967,0.1826034150449235888667637,0.1691565193950025381893121,0.1495959888165767320815017,0.1246289712555338720524763,0.0951585116824927848099251,0.0622535239386478928628438,0.0271524594117540948517806};


    /* n = 8 */
    static double x[4] = {0.1834346424956498049394761,0.5255324099163289858177390,0.7966664774136267395915539,0.9602898564975362316835609};
    static double w[4] = {0.3626837833783619829651504,0.3137066458778872873379622,0.2223810344533744705443560,0.1012285362903762591525314};

    double A,B,C,D,Ax,Cy,s;
    int i,j, m;

    m = (n+1)>>1;

    A = 0.5*(b-a);
    B = 0.5*(b+a);
    C = 0.5*(d-c);
    D = 0.5*(d+c);

    s = 0.0;
    for (i=0;i<m;i++)
    {
        Ax = A*x[i];
        for (j=0;j<m;j++)
        {
            Cy = C*x[j];
            s += w[i]*w[j]*( pointSpreadFunction(B+Ax,D+Cy,data)+pointSpreadFunction(Ax+B,D-Cy,data)
                             +pointSpreadFunction(B-Ax,D+Cy,data)+pointSpreadFunction(B-Ax,D-Cy,data));
        }
    }

    return C*A*s;
}
