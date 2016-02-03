#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>

//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
void step (cmplx* const psi, const double omega, const double dt, const double dx, const double xmin, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
  const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax - xmin)/(Nx-1);
	const double dt = 0.1*dx;
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

  const double omega = 0.2;
	const double lambda = 10;//2*M_PI/(omega*omega);

  const double alpha = sqrt(omega);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
			step (psi0, omega, dt, dx, xmin, Nx);


         	t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

	return 0;
}
//-----------------------------------
void step (cmplx* const psi, const double omega, const double dt, const double dx, const double xmin, const int Nx) {

	cmplx* psi_temp = new cmplx[Nx];
	cmplx* d = new cmplx[Nx];

	cmplx* o = new cmplx[Nx];
	cmplx* r = new cmplx[Nx];

	// entry for off-diagonal entries
	const cmplx a = cmplx(0.0, -dt/(4*dx*dx));

	// new vector for diagonal entries, fill after initializing
	for (int i=0; i<Nx; i++) {
		d[i] = cmplx(1.0, dt/(2*dx*dx) + dt*omega*omega/4.0*pow(xmin+i*dx, 2));
	}

	// fill temporary vector
	psi_temp[0]    = conj(d[0])*psi[0] + conj(a)*psi[1];
	psi_temp[Nx-1] = conj(a)*psi[Nx-2] + conj(d[Nx-1])*psi[Nx-1];

	for (int i=1; i<Nx-1; i++) {
		psi_temp[i] = conj(a)*psi[i-1] + conj(d[i])*psi[i] + conj(a)*psi[i+1];
	}

	// solve tridiagonal system
	o[0] = a/d[0];
	r[0] = psi_temp[0]/d[0];

	for (int i=1; i<Nx; i++) {
		o[i] = a/(d[i] - a*o[i-1]);
		r[i] = (psi_temp[i] - a*r[i-1])/(d[i] - a*o[i-1]);
	}

	// back substitution
	psi[Nx-1] = r[Nx-1];

	for (int i=2; i<=Nx; i++) {
		psi[Nx-i] = r[Nx-i] - o[Nx-i]*psi[Nx-i+1];
	}

	delete[] psi_temp;
	delete[] o;
	delete[] d;
	delete[] r;
}

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2.0 );
	}
}
