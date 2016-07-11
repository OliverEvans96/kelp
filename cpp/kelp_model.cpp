// kelp_model.cpp
// Oliver Evans
// Clarkson University REU 2016
// Created: Wed 06 Jul 2016 01:11:25 PM EDT
// Last Edited: Mon 11 Jul 2016 12:41:27 AM EDT

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

///////////////
// CONSTANTS //
///////////////

// Pi
const double PI = 3.141592653589793;

// Frond shape & ratio parameters
const double FS = 0.5;
const double FR = 2.0;
const double ALPHA = atan((1+FS)/(2*FR*FS));

// Vertical growth density (fronds/meter)
const double RHO = 10;


////////////////
// PROTOTYPES //
////////////////

// MATH FUNCTIONS //
double sign(double xx);
double min(double xx, double yy);
double max(double xx, double yy);
double factorial(int nn);
double I0(double xx,int nn);
double I0(double xx);
inline double vonMisesPDF(double xx,double kk,double mu);

// POLAR COORDINATE FUNCTIONS //
double theta_xy(double xx,double yy);
void polar_shift(double theta,double rr,double dz,double phi_s,double theta_s,double &theta_hat,double &rr_hat);

// NUMERICAL METHODS //
double csimp(double (*ff)(double),double aa,double bb,int nn);
double csimp(double (*ff)(double,double),double aa,double bb,int nn,double pp);
double csimp(double (*ff)(double,double,double),double aa,double bb,int nn,double pp1,double pp2);

// INTERPOLATION //
double linear(double* xx,double* yy,int nn,double xval);
double bilinear(double* xx,double* yy,double** zz,int nn,int mm,double xval,double yval);

// UTILITY FUNCTIONS //
inline int findBin(double aa, double bb,int nn,double xx);
inline int upperBin(double aa, double bb,int nn,double xx);
void print3d(double*** v, int nx, int ny,int nz, const char *name);
void write3d(double*** v, int nx, int ny,int nz, const char *name,ofstream &outFile);
void print2d(double** v, int nx, int ny, const char *name);
void print1d(double* v, int nn, const char *name);
void contourf(double* theta,double* rr,double** zz,int ntheta,int nr);

// MODEL-SPECIFIC FUNCTIONS //
inline double P_theta_f(double xx,double v_w,double theta_w);
inline double L_min_shade(double theta_p,double r_p,double theta_f);
inline double theta_min_shade(double theta_p,double r_p,double LL);
inline double theta_max_shade(double theta_p,double r_p,double LL);
double E_shade_3d(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L);
double P_shade_2d(double theta_p,double r_p,double v_w,double theta_w,int nLBin,double* LBin_vals,double* P_L);
double availableLight(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double surfaceIntensity,double attenuationCoef,double absorptionCoef);
inline void frondTransform(double ss,double tt,double &xx,double &yy,double theta_f,double LL);
inline double JJ(double ss,double tt,double LL);
double calculateLightAbsorption(double theta_f,double LL,double z_f,double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double surfaceIntensity,double attenuationCoef,double absorptionCoef);
inline double newLength(double lightAbsorbed,double LL,double tt);
void recalculateLengthDistribution(double v_w,double z_f,double nLBin,double dLBin,double LBin_min,double LBin_max,double LBin_vals,double* P_L,double tt,int nLBins);


//////////
// MAIN //
//////////

int main()
{
	// Bin convention:
	// Values stored in *Bin_vals are closed lower edges of bins
	// Values stored in lengthDist[ii][jj] are the proportion of the
	// population of fronds whose depth is in the interval
	// [ zBin_vals[ii] , zBin_vals[ii] )
	// whose length is in the interval
	// [ LBin_vals[ii] , LBin_vals[ii] )
	// (Upper half open intervals)

	// n*Bin is the number of bins (not number of edges)

	// Output to python
	ofstream outFile("../python/kelp_output.py");
	outFile << "from numpy import array" << endl;
	outFile << "import pickle" << endl;

	// Partition depth (units are meters)
	double dzBin = 0.5;
	double zBin_min = 0;
	double zBin_max = 10;
	int nzBin = int(floor((zBin_max-zBin_min)/dzBin)); // Length of zBin_vals
	double* zBin_vals = new double[nzBin];
	for(int ii=0;ii<nzBin;ii++)
		zBin_vals[ii] = zBin_min + dzBin*ii;
	
	// Partition length (units are meters)
	// First bin [-dLBin,0] is for dead fronds
	double dLBin = 0.1;
	double LBin_min = -dLBin;
	double LBin_max = 5;
	int nLBin = int(floor((LBin_max-LBin_min)/dLBin));
	double* LBin_vals = new double[nLBin];
	for(int ii=0;ii<nLBin;ii++)
		LBin_vals[ii] = LBin_min + dLBin*ii;

	// Population distribution values
	// Each depth is treated as a separate population
	// Each row should sum to 1
	double** lengthDist = new double*[nzBin];
	double sum;
	for(int ii=0;ii<nzBin;ii++)
	{
		lengthDist[ii] = new double[nLBin];

		// Initialize population to all be as small as possible
		// (But not dead)
		for(int jj = 0;jj<nLBin;jj++)
			lengthDist[ii][jj] = 0;
		lengthDist[ii][0] = 1;

		// Assign values
		sum = 0;
		for(int jj=1;jj<nLBin;jj++)
		{
			lengthDist[ii][jj] = (nzBin-ii)*exp(-pow(LBin_vals[jj]-3,2));
			//lengthDist[ii][jj] = 1;
			sum += lengthDist[ii][jj];
		}

		//Normalize
		for(int jj=0;jj<nLBin;jj++)
		{
			lengthDist[ii][jj] /= sum;
		}
	}

	/*
	// Only fronds in top layer
	lengthDist[0][0] = 0;
	lengthDist[0][nLBin-1] = 1;

	// And a layer in the middle
	lengthDist[5][0] = 0;
	lengthDist[5][nLBin-1] = 1;
	*/

	// Set up xy grid
	double xmin = -5;
	double xmax = 5;
	double dx = 1;
	int nx = int(floor((xmax-xmin)/dx));
	double* xx = new double[nx];
	for(int ii=0;ii<nx;ii++)
		xx[ii] = xmin + ii*dx;

	double ymin = -5;
	double ymax = 5;
	double dy = 1;
	int ny = int(floor((ymax-ymin)/dy));
	double* yy = new double[ny];
	for(int jj=0;jj<ny;jj++)
		yy[jj] = ymin + jj*dx;

	// Allocate 3d arrays
	double*** xx_3d = new double**[nx];
	double*** yy_3d = new double**[nx];
	double*** zz_3d = new double**[nx];
	double*** PP_3d = new double**[nx];
	for(int ii=0;ii<nx;ii++)
	{
		xx_3d[ii] = new double*[ny];
		yy_3d[ii] = new double*[ny];
		zz_3d[ii] = new double*[ny];
		PP_3d[ii] = new double*[ny];
		for(int jj=0;jj<ny;jj++)
		{
			xx_3d[ii][jj] = new double[nzBin];
			yy_3d[ii][jj] = new double[nzBin];
			zz_3d[ii][jj] = new double[nzBin];
			PP_3d[ii][jj] = new double[nzBin];
		}
	}

	double theta_p;
	double r_p;
	double theta_w = PI/2;
	double v_w = 1.0;
	double phi_s = 0;
	double theta_s = -PI/4;

	cout << "nzBin = " << nzBin << endl;
	//cout << "A = " << A << endl;

	// Calculate expected number of plants shading at each point in 3d grid
	for(int ii=0;ii<nx;ii++)
	{
		for(int jj=0;jj<ny;jj++)
		{
			theta_p = theta_xy(xx[ii],yy[jj]);
			r_p = sqrt(xx[ii]*xx[ii] + yy[jj]*yy[jj]);
			for(int kk=0;kk<nzBin;kk++)
			{
				//cout << "(" << ii << "," << jj << "," << kk << ")" << endl;
				xx_3d[ii][jj][kk] = xx[ii];
				yy_3d[ii][jj][kk] = yy[jj];
				zz_3d[ii][jj][kk] = zBin_vals[kk];
				PP_3d[ii][jj][kk] = E_shade_3d(theta_p,r_p,zBin_vals[kk],v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,lengthDist);
			}
		}
	}

	theta_p = PI/2;
	r_p = 1;

	// Write variables to file
	write3d(xx_3d,nx,ny,nzBin,"xx_3d",outFile);
	write3d(yy_3d,nx,ny,nzBin,"yy_3d",outFile);
	write3d(zz_3d,nx,ny,nzBin,"zz_3d",outFile);
	write3d(PP_3d,nx,ny,nzBin,"PP_3d",outFile);

	// Close output file
	outFile.close();

	// Run newly created python script
	system("mkdir -p ../python/kelp_pickle");
	system("python ../python/kelp_output.py");

	// Delete arrays
	delete [] xx;
	delete [] yy;
	for(int ii=0;ii<nx;ii++)
	{
		for(int jj=0;jj<ny;jj++)
		{
			delete [] xx_3d[ii][jj];
			delete [] yy_3d[ii][jj];
			delete [] zz_3d[ii][jj];
			delete [] PP_3d[ii][jj];
		}
		delete [] xx_3d[ii];
		delete [] yy_3d[ii];
		delete [] zz_3d[ii];
		delete [] PP_3d[ii];
	}
	delete [] xx_3d;
	delete [] yy_3d;
	delete [] zz_3d;
	delete [] PP_3d;

	for(int ii=0;ii<nzBin;ii++)
		delete lengthDist[ii];
	delete [] lengthDist;
	delete [] zBin_vals;
	delete [] LBin_vals;

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

// Minimum of two numbers
double min(double xx, double yy)
{
	if(xx<yy)
		return xx;
	else
		return yy;
}

// Maximum of two numbers
double max(double xx, double yy)
{
	if(xx>yy)
		return xx;
	else
		return yy;
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
// xx: x value
// nn: number of terms
// return: I0(xx)
double I0(double xx,int nn)
{
	double result = 0;
	
	for(int kk=0;kk<nn;kk++)
		result += pow((xx*xx/4),kk)/(pow(factorial(kk),2));

	return result;
}
// 20 terms seems like a reasonable default
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
inline double vonMisesPDF(double xx,double kk,double mu)
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

// Shift polar coordinate values along a cartesian line *towards* the sun
// theta: initial angular position
// rr: initial radial position
// dz: vertical distance *upward* to shift
// phi_s: zenith angle of the sun
// theta_s: azimuth angle of the sun
// theta_hat: new angular position
// rr_hat: new radial position
// * NOT YET CONSIDERING CHANGE IN PHI_S DUE TO REFRACTION *
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
	rr_hat = sqrt(xx_hat*xx_hat + yy_hat*yy_hat);
}


///////////////////////
// NUMERICAL METHODS //
///////////////////////

// Composite Simpson's rule for integration (1d)
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
// Allow for two double parameters, pp1 and pp2, to be passed to ff
double csimp(double (*ff)(double,double,double),double aa,double bb,int nn,double pp1,double pp2)
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
		total += hh/6 * (ff(x_i,pp1,pp2) + 4*ff((x_i+x_ip1)/2,pp1,pp2) + ff(x_ip1,pp1,pp2));
	}

	return total;
}


///////////////////
// INTERPOLATION //
///////////////////

// Linear interpolation
// xx: 1d ordered arrays of known x values
// yy: 1d ordered array of known y values
// nn: length of xx and yy
// xval: x value to interpolate
// return: interpolated y value
double linear(double* xx,double* yy,int nn,double xval)
{
	int ii = findBin(xx[0],xx[nn-1],nn,xval);
	if(ii < nn)
		return yy[ii] + (xval-xx[ii])/(xx[ii+1]-xx[ii]) * (yy[ii+1] - yy[ii]);
	else
		return yy[nn-1];
}

// Bilinear interpolation
// xx: 1d ordered array of known x values
// yy: 1d ordered array of known y values
// zz: 2d array of known z values
// nn: length of xx
// mm: length of yy
// xval: x coordinate of interpolation point
// yval: y coordinate of interpolation point
// return: interpolated z value
double bilinear(double* xx,double* yy,double** zz,int nn,int mm,double xval,double yval)
{
	// Find closest known points
	int ii = findBin(xx[0],xx[nn-1],nn,xval);
	int jj = findBin(yy[0],yy[mm-1],mm,yval);

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
		cout << "Out of bounds!" << endl;
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
inline int findBin(double aa, double bb,int nn,double xx)
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
inline int upperBin(double aa, double bb,int nn,double xx)
{
	return ceil((xx-aa)/(bb-aa)*nn);
}

// Print 3d array for input to python
void print3d(double*** v, int nx, int ny,int nz, const char *name)
{
	cout << name << " = array([";
	for(int ii=0;ii<nx;ii++) 
	{
		cout << "[";
		for(int jj=0;jj<ny;jj++)
		{
			cout << "[";
			for(int kk=0;kk<nz;kk++)
			{
				cout << v[ii][jj][kk];
				if(kk<nz-1)
					cout << ",";
			}
			cout << "]";
			if(jj<ny-1)
				cout << "," << endl;
		}
		cout << "]";
		if(ii<nx-1)
			cout << "," << endl << endl;
	}
	cout << "])" << endl;
}

// Write 3d array for input to python to file
void write3d(double*** v, int nx, int ny,int nz, const char *name,ofstream &outFile)
{
	outFile << name << " = array([";
	for(int ii=0;ii<nx;ii++) 
	{
		outFile << "[";
		for(int jj=0;jj<ny;jj++)
		{
			outFile << "[";
			for(int kk=0;kk<nz;kk++)
			{
				outFile << v[ii][jj][kk];
				if(kk<nz-1)
					outFile << ",";
			}
			outFile << "]";
			if(jj<ny-1)
				outFile << "," << endl;
		}
		outFile << "]";
		if(ii<nx-1)
			outFile << "," << endl << endl;
	}
	outFile << "])" << endl;

	// Save with pickle
	outFile << "with open('../python/kelp_pickle/" << name << ".pickle','wb') as pickle_file:" << endl;
	outFile << "    pickle.dump(" << name << ",pickle_file)" << endl;
}



// Print 2d array for input to python
void print2d(double** v, int nx, int ny, const char *name)
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

// Print 1d array for input to python
void print1d(double* v, int nn, const char *name)
{
	cout << name << " = array([";
	for(int ii=0;ii<nn-1;ii++)
		cout << v[ii] << ",";
	cout << v[nn-1] << "])" << endl;
}

// Create polar filled contour plot of values in Python
void contourf(double* theta,double* rr,double** zz,int ntheta,int nr)
{
	print1d(theta,ntheta,"th");
	print1d(rr,nr,"rr");
	print2d(zz,ntheta,nr,"zz");
	cout << "rr,th = meshgrid(rr,th)" << endl;
	cout << "clf()" << endl;
	cout << "gca(projection='polar')" << endl;
	cout << "contourf(th,rr,zz,cmap='YlGn')" << endl;
	cout << "colorbar()" << endl;
}


//////////////////////////////
// MODEL-SPECIFIC FUNCTIONS //
//////////////////////////////

// PDF of frond angle distribution
// v_w: water current speed
// theta_w: angle of water current in radians
// xx: input value
// return: probability density at xx
inline double P_theta_f(double xx,double v_w,double theta_w)
{
	double result = 1/(2*PI);
	if(v_w != 0)
		result = vonMisesPDF(xx,v_w,theta_w);
	return result;
}

// Minimum LL (frond length) for shading a given point
// as a function of theta_f (frond angle)
inline double L_min_shade(double theta_p,double r_p,double theta_f)
{
	double theta_prime = theta_p - theta_f + PI/2;
	double S = sign(PI/2-theta_prime);
	return r_p*(sin(theta_prime)+2*S*FR/(1+FS)*cos(theta_prime));
}

// Minimum theta_f (frond angle ) for shading a given point
// as a function of LL (frond length)
// theta_max_shade = 2*theta_p - theta_min_shade
inline double theta_min_shade(double theta_p,double r_p,double LL)
{
	double qq = 2*FR/(1+FS);
	// Unshifted value, assuming point is at theta_p = PI/2
	double unshifted = min( atan(1/qq) + acos(LL/r_p/sqrt(1+qq*qq)) , PI/2 );
	// Rotate point to theta_p
	return unshifted + theta_p - PI/2;
}

// Expected number of fronds shading a point in 3d space
// considering only incident light (not diffuse)
// Integrate P_shade_2d over z, shifting theta_p and r_p towards the sun appropriately
double E_shade_3d(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L)
{
	// Probability of shading
	double PP = 0;

	// Shifted polar coordinates
	double theta_hat,r_hat;

	// Current depth
	double z_f;

	// Index of z bin containing point
	int k_p = findBin(zBin_min,zBin_max,nzBin,z_p);

	// Loop over depths from surface to z_p, inclusive
	for(int kk=0;kk<k_p;kk++)
	{
		z_f = zBin_vals[kk];
		polar_shift(theta_p,r_p,z_p-z_f,phi_s,theta_s,theta_hat,r_hat);
		PP += RHO/dzBin * P_shade_2d(theta_hat,r_hat,v_w,theta_w,nLBin,LBin_vals,P_L[kk]);
	}

	/*
	print1d(z_f,nz,"z_f");
	print1d(PP_z,nz,"PP_z");
	cout << "plot(z_f,PP_z)" << endl;
	*/

	return PP;
}

double P_shade_2d(double theta_p,double r_p,double v_w,double theta_w,int nLBin,double* LBin_vals,double* P_L)
{
	// Recalculate information about dLBin_vals
	double dLBin = LBin_vals[1] - LBin_vals[0];
	double LBin_min = LBin_vals[0];
	double LBin_max = LBin_vals[nLBin-1] + dLBin;

	// Whether integration has to stop in theta-variable region
	bool early_cutoff = false;

	// Calculate LL limits of numerical integration
	double L_star = L_min_shade(theta_p,r_p,theta_p+ALPHA);
	int k0 = findBin(LBin_min,LBin_max,nLBin,r_p);
	// Don't let k_star exceed number of bins
	int k_star = upperBin(LBin_min,LBin_max,nLBin,L_star);

	// Check for early cutoff
	if(nLBin <= k_star)
	{
		early_cutoff = true;
		k_star = nLBin;
	}

	// theta_f limits of integration for this length & next & avg
	double theta_min_k;
	double theta_min_kp1;
	double theta_min,theta_max;

	// Number of sub-intervals to use for integrating P_theta_f
	int n_csimp = 20;

	// Probability of shading in theta-variable, theta-constant regions and total
	double P_var = 0;
	double P_const = 0;
	double P_total = 0;

	// Integrate very first bin separately
	// Don't let theta_min exceed theta_p
	// Left endpoint
	theta_min_k = theta_min_shade(theta_p,r_p,LBin_vals[k0]);
	// Right endpoint
	theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[k0+1]);
	// Use average of left and right endpoints for integration limits
	theta_min = (theta_min_k + theta_min_kp1)/2;
	theta_max = 2*theta_p - theta_min;

	// Accumulate integral
	P_var += csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w) * P_L[k0];

	// Integrate over remainder of theta-variable region
	for(int kk=k0+1;kk<k_star-1;kk++)
	{
		// Left endpoint
		theta_min_k = theta_min_kp1;

		// Right endpoint
		theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[kk+1]);

		// Use average of left and right endpoints for integration limits
		theta_min = (theta_min_k + theta_min_kp1)/2;
		theta_max = 2*theta_p - theta_min;

		// Accumulate integral
		P_var += csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w) * P_L[kk];
	}

	// Integrate over last theta-variable bin separately to efficiently check for early cutoff
	// Left endpoint
	theta_min_k = theta_min_kp1;

	// Right endpoint
	if(early_cutoff)
		theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[k_star-1]+dLBin);
	else
		theta_min_kp1 = theta_min_shade(theta_p,r_p,LBin_vals[k_star]);

	// Use average of left and right endpoints for integration limits
	theta_min = (theta_min_k + theta_min_kp1)/2;
	theta_max = 2*theta_p - theta_min;

	// Accumulate integral
	P_var += csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w) * P_L[k_star-1];

	// Sum lengths
	theta_min = theta_p - ALPHA;
	theta_max = theta_p + ALPHA;
	for(int kk=k_star;kk<nLBin;kk++)
		P_const += P_L[kk];

	// Multiply by integral over theta_f
	P_const *= csimp(P_theta_f,theta_min,theta_max,n_csimp,v_w,theta_w);

	// Calculate total integral
	P_total = P_var + P_const;

	return P_total;
}

// Calculate available light at a given 3d point
// Not considering ambient light
double availableLight(double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double surfaceIntensity,double attenuationCoef,double absorptionCoef)
{
	double nShade = E_shade_3d(theta_p,r_p,zBin_vals[kk],v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,lengthDist);
	double absorptionFactor = pow((1-absorbtion_coef),nShade);
	double attenuationFactor = exp(-attenuationCoef*z_p);

	return surfaceIntensity * absorptionFactor * attenuationFactor;
}

// Transform a point (ss,tt) in the unit square: [-1,1]x[-1,1] to the corresponding point (xx,yy) on a frond of length L at an angle theta_f
inline void frondTransform(double ss,double tt,double &xx,double &yy,double theta_f,double LL)
{
	xx = LL*(-FR*(FS*(s - 1)*(t + 1) - (s + 1)*(2*FS + t + 1))*cos(theta_f) 
			+ (FS + 1)*(s - t)*sin(theta_f))/(4*FR*(FS + 1)) 
	yy = -LL*(FR*(FS*(s - 1)*(t + 1) - (s + 1)*(2*FS + t + 1))*sin(theta_f) 
			+ (FS + 1)*(s - t)*cos(theta_f))/(4*FR*(FS + 1))
}

// Jacobian of frond transform
inline double JJ(double ss,double tt,double LL)
{
	return L**2*(-FS**2*ss - FS**2*tt + 2*FS**2 + 4*FS + ss + tt + 2)/(16*FR*(FS + 1)**2);
}

// Calculate the light absorbed by a frond with a particular angle (theta_f) and length (LL)
// at a particular depth (z_f)
// Use n=2 Gaussian quadrature product rule to integrate light field over frond area by transforming frond to unit square
double calculateLightAbsorption(double theta_f,double LL,double z_f,double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double surfaceIntensity,double attenuationCoef,double absorptionCoef)
{
	// Abscissas on unit square
	static double sa[2] = {-1/sqrt(3),1/sqrt(3)};
	
	// Weights (not necessary for n=2)
	// static double ww[2] = {1,1};

	// Coordinates of frond abscissae
	double xfa,yfa,tfa,rfa;

	// Total light absorbed by frond
	double lightAbsorbed = 0;

	// Part of integral due to a specific abscissa
	double absc_val;

	// Loop through 2x2 abscissas
	for(int ii=0;ii<2;ii++)
	{
		for(int jj=0;jj<2;jj++)
		{
			// Calculate abscissae
			frondTransform(sa[ii],sa[jj],xfa,yfa,theta_f,LL);

			// Convert to polar coordinates
			tfa = theta_xy(xfa,xya);
			rfa = sqrt(xfa*xfa + yfa*yfa);

			// Calculate contribution of this abscissa
			absc_val = availableLight(theta_p,r_p,z_p,v_w,theta_w,phi_s,theta_s,gnLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L,surfaceIntensity,attenuationCoef,absorptionCoef);

			// Multiply by weights according to quadrature product rule
			// absc_val *= ww[ii] * ww[jj];

			// Multiply by the Jacobian and add to total
			lightAbsorbed += absc_val * abs(JJ(sa[ii],sa[jj],LL));
		}
	}

}

// New frond length from growth as a function of light absorbed, frond length, and time
inline double newLength(double lightAbsorbed,double LL,double tt)
{
	double light_param = 1;
	double length_param = 1;
	double time_param = 1;
	double growth = exp(lightAbsorbed*light_param - LL*length_param - tt*time_param);
	return LL + growth;
}

// Recalculate length distribution by integrating P_theta_f*P_L over R_i for each L_i,
// where R_i is the region such that newLength(theta_f,LL) \in [ L_i , L_{i+1} )
void recalculateLengthDistribution(double v_w,double z_f,double nLBin,double dLBin,double LBin_min,double LBin_max,double LBin_vals,double* P_L,double tt,int nLBins)
calculateLightAbsorption(double theta_f,double LL,double z_f,double theta_p,double r_p,double z_p,double v_w,double theta_w,double phi_s,double theta_s,int nLBin,double* LBin_vals,int nzBin,double dzBin,double zBin_min,double zBin_max,double* zBin_vals,double** P_L,double surfaceIntensity,double attenuationCoef,double absorptionCoef)
{
	// Number of points to sample from theta_f
	int n_theta_f = 20;
	double dtheta_f = 2*PI/n_theta_f_points;

	// Current sample coordinates
	double theta_f = -PI;
	double LL_current;

	// Light absorbed by frond
	double lightAbsorbed;

	// theta_f and LL coordinates of points in R_i for a given LL bin
	vector<double>* R_i_theta_f = new vector<double>[nLBin];
	vector<double>* R_i_LL = new vector<double>[nLBin];

	// Number of points which fall into each LL bin
	vector<int>* nPoints = new vector<int>[nLBin];

	// New length coordinate
	double LL_new;

	// Bin of new length coordinate
	int new_bin;

	// New length distribution
	double* P_L_new = new double[nLBin];

	// Loop through theta_f-LL space
	for(int ii=0;ii<n_theta_f;ii++)
	{
		for(int jj=0;jj<nLBin;jj++)
		{
			// Get LL value for this index (use bin midpoint)
			double LL_current = LBin_vals[jj] + dLBin/2;

			// Calculate light absorbed by a frond with this particular theta_f and LL
			lightAbsorbed = calculateLightAbsorption(theta_f,LL,z_f,theta_p,r_p,z_p,v_w,theta_w,phi_s,theta_s,nLBin,LBin_vals,nzBin,dzBin,zBin_min,zBin_max,zBin_vals,P_L,surfaceIntensity,attenuationCoef,absorptionCoef)

			// Calculate new frond size for this type of frond
			LL_new = newLength(lightAbsorbed,LL_current,tt);

			// Determine which bin this new length falls into
			new_bin = findBin(LBin_min,LBin_max,nLBin,LL_new);

			// Save the indices of this theta_f,LL pair in the appropriate vector
			R_i_theta_f[new_bin].push_back(LL_current);
			R_i_LL[new_bin].push_back(theta_f);

			// Increment point counter for this bin
			nPoints[new_bin]++;
		}

		// Increment sample point
		theta_f += dtheta_f;
	}

	// Integrate 2d population distribution over appropriate region for each new length bin
	// Loop through new length bins
	for(int kk=0;kk<nLBin;kk++)
	{
		// Integrate 2d population distribution over theta_f,LL points which map to this bin
		// This is the (approximate) percentage of the population 
		// which will fall in this bin after growth (next time step)

		// Loop over points (use 2d midpoint rule)
		for(int nn=0;nn<nPoints[kk];nn++)
			P_L_new[kk] += P_theta_f(R_i_theta_f[kk][nn],v_w,theta_w)*dtheta_f*P_L(R_i_LL[kk][nn]);
	}

	// Update length distribution
	for(int kk=0;kk<nLBin;kk++)
		P_L[kk] = P_L_new[kk];
}


/////////
// FIN //
/////////

