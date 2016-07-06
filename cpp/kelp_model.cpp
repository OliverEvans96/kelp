// kelp_model.cpp
// Oliver Evans
// Clarkson University REU 2016
// Created: Wed 06 Jul 2016 01:11:25 PM EDT
// Last Edited: Wed 06 Jul 2016 07:27:35 PM EDT

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double PI = 3.141592653589793;

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
		//x_{i+1}
		x_ip1 = aa + hh*(ii+1);
		total += hh/6 * (ff(x_i) + 4*ff((x_i+x_ip1)/2) + ff(x_ip1));
	}

	return total;
}

// Factorial function
double factorial(int nn)
{
	double result = 1;
	for(int ii=1;ii<=nn;ii++)
		result *= ii;
	return result;
}

// Bin finder
// Given an interval between aa and bb divided into nn bins,
// find the number of the bin which xx belongs to.
// aa: min value
// bb: max value
// nn: number of bins
// return: index of bin that xx belongs to
int findBin(double aa, double bb,int nn,double xx)
{
	return floor((xx-aa)/(bb-aa)*nn);
}

// Print 2D array for input to python
void print2d(double **v, int nx, int ny, const char *name)
{
	cout << name << " = array([";
	for(int ii=0;ii<nx-1;ii++) 
	{
		cout << "[";
		for(int jj=0;jj<ny-1;jj++) 
			cout << v[ii][jj] << ",";
		cout << v[ii][ny-1] << "],";
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

double ff(double xx,double yy)
{
	return sin(PI*xx)+sin(PI*yy);
}

int main()
{
	// Partition depth (units are meters)
	double zmin = 0;
	double zmax = 1;
	double dz = 0.5;
	int nz = int(floor((zmax-zmin)/dz))+1; // Length of zvals
	double *zvals = new double[nz];
	for(int ii=0;ii<nz;ii++)
		zvals[ii] = zmin + dz*ii;
	
	// Partition length (units are meters)
	double Lmin = 0;
	double Lmax = 2;
	double dL = 0.5;
	int nL = int(floor((Lmax-Lmin)/dL))+1;
	double *Lvals = new double[nL];
	for(int ii=0;ii<nL;ii++)
		Lvals[ii] = Lmin + dL*ii;

	// Test values
	double **Pvals = new double*[nz];
	for(int ii=0;ii<nz;ii++)
	{
		Pvals[ii] = new double[nL];
		for(int jj = 0;jj<nL;jj++)
			Pvals[ii][jj] = ff(zvals[ii],Lvals[jj]);
	}

	double dxx = 0.1;
	int nxx = int(floor((zmax-zmin)/dxx))+1;
	double *xx = new double[nxx];
	for(int ii=0;ii<nxx;ii++)
		xx[ii] = zmin + dxx*ii;

	double dyy = 0.1;
	int nyy = int(floor((Lmax-Lmin)/dyy))+1;
	double *yy = new double[nyy];
	for(int ii=0;ii<nyy;ii++)
		yy[ii] = Lmin + dyy*ii;

	// Interpolated values
	double **Interpvals = new double*[nxx];
	for(int ii=0;ii<nxx;ii++)
	{
		Interpvals[ii] = new double[nyy];
		for(int jj = 0;jj<nyy;jj++)
			Interpvals[ii][jj] = bilinear(zvals,Lvals,Pvals,nz,nL,xx[ii],yy[jj]);
	}

	print1d(zvals,nz,"z_vals");
	print1d(Lvals,nL,"L_vals");
	print2d(Pvals,nz,nL,"p_vals");

	print1d(xx,nxx,"xx");
	print1d(yy,nyy,"yy");
	print2d(Interpvals,nxx,nyy,"interp");

	// Delete arrays
	for(int ii=0;ii<nz;ii++)
		delete Pvals[ii];
	delete [] Pvals;
	delete zvals;
	delete Lvals;

	return 0;
}
