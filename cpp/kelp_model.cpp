// kelp_model.cpp
// Oliver Evans
// Clarkson University REU 2016
// Created: Wed 06 Jul 2016 01:11:25 PM EDT
// Last Edited: Thu 07 Jul 2016 11:26:16 AM EDT

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

///////////////
// CONSTANTS //
///////////////

// Pi
const double PI = 3.141592653589793;

// Frond shape & ratio parameters
const double FS = 0.5;
const double FR = 2.0;


////////////////
// PROTOTYPES //
////////////////

// MATH FUNCTIONS //
double sign(double xx);
double factorial(int nn);
double I0(double xx);
double vonMisesPDF(double kk,double mu,double xx);

// POLAR COORDINATE FUNCTIONS //
double theta_xy(double xx,double yy);
void polar_shift(double theta,double rr,double dz,double phi_s,double theta_s,double &theta_hat,double &rr_hat);

// NUMERICAL INTEGRATION //
double csimp(double (*ff)(double),double aa,double bb,int nn);
double csimp(double (*ff)(double,double),double aa,double bb,int nn,double pp);

// INTERPOLATION //
double bilinear(double *xx,double *yy,double **zz,int nn,int mm,double xval,double yval);

// UTILITY FUNCTIONS //
int findBin(double aa, double bb,int nn,double xx);
int upperBin(double aa, double bb,int nn,double xx);
void print2d(double **v, int nx, int ny, const char *name);
void print1d(double *v, int nn, const char *name);

// MODEL-SPECIFIC FUNCTIONS //
double P_theta_f(double vw,double theta_w_deg,double xx);
double Lmin(double theta_p,double r_p,double theta_f);
double P_shade(double theta_p,double r_p,double z_p);
double P_shade_2D(double theta_p,double r_p,int nL,double dL);


//////////
// MAIN //
//////////

int main()
{
	// Calculate frond angle parameter
	double alpha = atan((1+FS)/(2*FR*FS));

	// Bin convention:
	// Values stored in *Bin_vals are closed lower edges of bins
	// Values stored in lengthDist[ii][jj] are the proportion of the
	// population of fronds whose depth is in the interval
	// [ zBin_vals[ii] , zBin_vals[ii] )
	// whose length is in the interval
	// [ LBin_vals[ii] , LBin_vals[ii] )
	// (Upper half open intervals)

	// n*Bin is the number of bins (not number of edges)

	// Partition depth (units are meters)
	double dzBin = 1;
	double zBin_min = 0;
	double zBin_max = 10;
	int nzBin = int(floor((zBin_max-zBin_min)/dzBin)); // Length of zBin_vals
	double *zBin_vals = new double[nzBin];
	for(int ii=0;ii<nzBin;ii++)
		zBin_vals[ii] = zBin_min + dzBin*ii;
	
	// Partition length (units are meters)
	// First bin [-dLBin,0] is for dead fronds
	double dLBin = 1;
	double LBin_min = -dLBin;
	double LBin_max = 5;
	int nLBin = int(floor((LBin_max-LBin_min)/dLBin));
	double *LBin_vals = new double[nLBin];
	for(int ii=0;ii<nLBin;ii++)
		LBin_vals[ii] = LBin_min + dLBin*ii;

	// Population distribution values
	// Each depth is treated as a separate population
	// Each row should sum to 1
	double **lengthDist = new double*[nzBin];
	for(int ii=0;ii<nzBin;ii++)
	{
		lengthDist[ii] = new double[nLBin];
		// Initialize population to all be as small as possible
		// (But not dead)
		lengthDist[ii][1] = 1;
		for(int jj = 2;jj<nLBin;jj++)
			lengthDist[ii][jj] = 0;
	}





	double test_val = -1;
	cout << "LBin_min = " << LBin_min << endl;
	cout << "LBin_max = " << LBin_max << endl;
	cout << "dLBin = " << dLBin << endl;
	cout << "nLBin = " << nLBin << endl;
	print1d(LBin_vals,nLBin,"LBin_vals");
	cout << endl;
	for(int ii=0;ii<10;ii++)
	{
		test_val += 1.0/3;
		cout << "Bin #s of (" << setw(9) << test_val << ") are: ";
		cout << findBin(LBin_min,LBin_max,nLBin,test_val) << " and ";
		cout << upperBin(LBin_min,LBin_max,nLBin,test_val) << endl;
	}

	cout << endl;
	cout << "zBin_min = " << zBin_min << endl;
	cout << "zBin_max = " << zBin_max << endl;
	cout << "dzBin = " << dzBin << endl;
	cout << "nzBin = " << nzBin << endl;
	print1d(zBin_vals,nzBin,"zBin_vals");









	// Delete arrays
	for(int ii=0;ii<nzBin;ii++)
		delete lengthDist[ii];
	delete [] lengthDist;
	delete zBin_vals;
	delete LBin_vals;

	return 0;
}

////////////////////
// MATH FUNCTIONS //
////////////////////

// Sign function
double sign(double xx)
{
	if(xx > 0)
		return 1;
	else if(xx < 0)
		return -1;
	else
		return 0;
}

// Factorial function
double factorial(int nn)
{
	double result = 1;
	for(int ii=1;ii<=nn;ii++)
		result *= ii;
	return result;
}

// Modified bessel function of the first kind of order 0
// Summation formula from Wolfram Mathworld: http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
double I0(double xx)
{
	double result = 0;
	
	// 20 terms seems reasonable
	for(int kk=0;kk<20;kk++)
		result += pow((xx*xx/4),kk)/(pow(factorial(kk),2));

	return result;
}

// von Mises Distribution PDF
// kk: kappa - sharpness parameter
// mu: horizontal shift
// xx: input value
// returns: PDF evaluated at xx
double vonMisesPDF(double kk,double mu,double xx)
{
	return exp(kk*cos(xx-mu))/(2*PI*I0(kk));
}


////////////////////////////////
// POLAR COORDINATE FUNCTIONS //
////////////////////////////////

// Calculate polar theta from x,y values
double theta_xy(double xx,double yy)
{
	double theta = sign(yy)*PI/2;
	if(xx!=0)
	{
		theta = atan(yy/xx);
		if(xx<0)
			theta += PI;
	}

	// Shift theta to the interval [-PI,PI]
	if(theta>=PI)
		theta -= 2*PI;

	return theta;
}

// Shift polar coordinate values along a cartesian line towards the sun
// theta: initial angular position
// rr: initial radial position
// dz: vertical distance to shift
// phi_s: zenith angle of the sun
// theta_s: azimuth angle of the sun
// theta_hat: new angular position
// rr_hat: new radial position
void polar_shift(double theta,double rr,double dz,double phi_s,double theta_s,double &theta_hat,double &rr_hat)
{
	// Initial cartesian coordinates
	double xx = rr*cos(theta);
	double yy = rr*sin(theta);

	// Cartesian shifts
	double dx = dz*tan(phi_s)*cos(theta_s);
	double dy = dz*tan(phi_s)*sin(theta_s);

	// New cartesian coordinates
	double xx_hat = xx + dx;
	double yy_hat = yy + dy;

	// New polar coordinates
	theta_hat = theta_xy(xx_hat,yy_hat);
	rr_hat = sqrt(xx*xx + yy*yy);
}


///////////////////////////
// NUMERICAL INTEGRATION //
///////////////////////////

// Composite Simpson's rule for integration (1D)
// ff: integrand
// aa: lower limit of integration
// bb: upper limit of integration
// nn: number of sub-intervals
// return: definite integral of f from a to b with n sub-intervals
double csimp(double (*ff)(double),double aa,double bb,int nn)
{
	double total = 0;
	double hh = (bb-aa)/nn;
	double x_i,x_ip1;

	// Loop through intervals
	for(int ii=0;ii<nn;ii++)
	{
		//x_i
		x_i = aa + hh*ii;
		//x{i+1}
		x_ip1 = aa + hh*(ii+1);
		total += hh/6 * (ff(x_i) + 4*ff((x_i+x_ip1)/2) + ff(x_ip1));
	}

	return total;
}
// Allow for one double parameter, pp, to be passed to ff
double csimp(double (*ff)(double,double),double aa,double bb,int nn,double pp)
{
	double total = 0;
	double hh = (bb-aa)/nn;
	double x_i,x_ip1;

	// Loop through intervals
	for(int ii=0;ii<nn;ii++)
	{
		//x_i
		x_i = aa + hh*ii;
		//x{i+1}
		x_ip1 = aa + hh*(ii+1);
		total += hh/6 * (ff(x_i,pp) + 4*ff((x_i+x_ip1)/2,pp) + ff(x_ip1,pp));
	}

	return total;
}


///////////////////
// INTERPOLATION //
///////////////////

// Bilinear interpolation
// xx: 1D ordered array of known x values
// yy: 1D ordered array of known y values
// zz: 2D array of known z values
// nn: length of xx
// mm: length of yy
// xval: x coordinate of interpolation point
// yval: y coordinate of interpolation point
// return: interpolated z value
double bilinear(double *xx,double *yy,double **zz,int nn,int mm,double xval,double yval)
{
	// Find closest known points
	int ii = findBin(xx[0],xx[nn-1],nn-1,xval);
	int jj = findBin(yy[0],yy[mm-1],mm-1,yval);

	double result;

	// Upper left
	double A = (xval-xx[ii])*(yy[jj+1]-yval);
	// Upper right
	double B = (xx[ii+1]-xval)*(yy[jj+1]-yval);
	// Lower left
	double C = (xval-xx[ii])*(yval-yy[jj]);
	// Lower right
	double D = (xx[ii+1]-xval)*(yval-yy[jj]);

	// Not on upper edges
	if(ii < nn-1 && jj < mm-1)
		result = (A*zz[ii+1][jj] + B*zz[ii][jj] + C*zz[ii+1][jj+1] + D*zz[ii][jj+1])/(A+B+C+D);

	// Right edge
	else if (ii == nn-1 && jj < mm-1)
	{
		A = yy[jj+1]-yval;
		C = yval-yy[jj];
		result = (A*zz[nn-1][jj] + C*zz[nn-1][jj+1])/(A+C);
	}

	// Upper edge
	else if (ii < nn-1 && jj == mm-1)
	{
		B = xx[ii+1]-xval;
		D = xx[ii+1]-xval;
		result = (B*zz[ii][mm-1] + D*zz[ii+1][mm-1])/(B+D);
	}

	// Upper left corner
	else if (ii == nn-1 && jj == mm-1)
		result = zz[nn-1][mm-1];

	// Out of bounds
	else
	{
		// cout << "Out of bounds!" << endl;
		result = 0;
	}
	return result;
}


///////////////////////
// UTILITY FUNCTIONS //
///////////////////////

// Bin finder
// Given an interval between aa and bb divided into nn bins,
// find the number of the bin which xx belongs to.
// aa: min value
// bb: max value
// nn: number of bins
// xx: given value
// return: index of bin that xx belongs to
int findBin(double aa, double bb,int nn,double xx)
{
	return floor((xx-aa)/(bb-aa)*nn);
}

// Index of the next bin if xx does not fall
// exactly on a bin edge, or the current bin if it does
// Remember: bins are upper half open.
// **This is not necessarily the bin that xx belongs to**
// aa: min value
// bb: max value
// nn: number of bins
// xx: given value
// return: index of bin
int upperBin(double aa, double bb,int nn,double xx)
{
	return ceil((xx-aa)/(bb-aa)*nn);
}

// Print 2D array for input to python
void print2d(double **v, int nx, int ny, const char *name)
{
	cout << name << " = array([" << endl;
	for(int ii=0;ii<nx-1;ii++) 
	{
		cout << "[";
		for(int jj=0;jj<ny-1;jj++) 
			cout << v[ii][jj] << ",";
		cout << v[ii][ny-1] << "]," << endl;
	}
	cout << "[";
	for(int jj=0;jj<ny-1;jj++) 
		cout << v[nx-1][jj] << ",";
	cout << v[nx-1][ny-1] << "]])" << endl;
}

// Print 1D array for input to python
void print1d(double *v, int nn, const char *name)
{
	cout << name << " = array([";
	for(int ii=0;ii<nn-1;ii++)
		cout << v[ii] << ",";
	cout << v[nn-1] << "])" << endl;
}


//////////////////////////////
// MODEL-SPECIFIC FUNCTIONS //
//////////////////////////////

// PDF of frond angle distribution
// vw: water current speed
// theta_w_deg: angle of water current in degrees
// xx: input value
// return: probability density at xx
double P_theta_f(double vw,double theta_w_deg,double xx)
{
	double theta_w = theta_w_deg * PI/180;
	double result = 1/(2*PI);
	if(vw != 0)
		result = vonMisesPDF(1/vw,theta_w,xx);
	return result;
}

// Minimum L value for shading a given point
double Lmin(double theta_p,double r_p,double theta_f)
{
	double theta_prime = theta_p - theta_f + PI/2;
	double S = sign(PI/2-theta_prime);
	return r_p*(sin(theta_prime)+2*S*FR/(1+FS)*cos(theta_prime));
}

// Probability of a polar point being shaded by kelp
// considering only incident light (not diffuse)
double P_shade(double theta_p,double r_p,double z_p)
{
	return 0;
}

double P_shade_2D(double theta_p,double r_p,int nL,double dL)
{
	int k0 = 0;
	int k_star = 1;

	double result = 0;

	// Integrate over theta-variable region
	for(int kk=0;kk<k_star;kk++)
	{
	}

	// Integrate over theta-constant region
	for(int kk=k_star;kk<nL;kk++)
	{
	}

	// Multiply by L-interval width
	result *= dL;
}


/////////
// FIN //
/////////

