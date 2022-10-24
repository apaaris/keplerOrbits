#include<iostream>
#include<cmath>
#include<omp.h>
#include<fstream>
#include<chrono>

# define E 0.9
# define NSTEPS 1e4
# define folder "ec9"
double T(){

	return std::pow((2.0 * M_PI) / (1.0 - E), (3.0 / 2.0));
}

double fx(double x, double y){

	double r = std::sqrt(x*x + y*y);
	return (-1) / std::pow(r,3) * x;

}

double fy(double x, double y){

	double r = std::sqrt(x*x + y*y);
	return (-1) / std::pow(r,3) * y;

}


void eEuler (double* x, double* y, double* vx, double* vy, double dt){


	for(size_t i = 0; i < NSTEPS; i++){
		x[i+1] = x[i] + vx[i] * dt;
		y[i+1] = y[i] + vy[i] * dt;

		vx[i+1] = vx[i] + fx(x[i],y[i]) * dt;
		vy[i+1] = vy[i] + fy(x[i],y[i]) * dt;
	}

	// FLOp:
	//
	// 28 N

}


void rk2 (double * x, double * y, double * vx, double * vy, double dt){

	double k1x, k1y;

	for(size_t i = 0; i < NSTEPS; i++){

		k1x = dt * fx(x[i], y[i]) / 2.0;
		k1y = dt * fy(x[i], y[i]) / 2.0;

		vx[i+1] = vx[i] + dt * fx(x[i] + k1x, y[i]);
		vy[i+1] = vy[i] + dt * fy(x[i], y[i] + k1y);

		x[i+1] = x[i] + vx[i] * dt;
		y[i+1] = y[i] + vy[i] * dt;

	}

	// FLOp:
	//
	// 54 N

}


void rk4 (double * x, double * y, double * vx, double * vy, double dt){

	double k1vx,k1rx,k2vx,k2rx,k3vx,k3rx,k4vx,k4rx;
	double k1vy,k1ry,k2vy,k2ry,k3vy,k3ry,k4vy,k4ry;

	for(size_t i = 0; i < NSTEPS; i++){
		k1vx = fx(x[i],y[i]);
		k1rx = vx[i];
		k2vx = fx(x[i] + (dt/2) * k1rx, y[i]);
		k2rx = vx[i] + (dt/2) * k1vx;
		k3vx = fx(x[i] + (dt/2) * k2rx, y[i]);
		k3rx = vx[i] + (dt/2) * k2vx;
		k4vx = fx(x[i] + dt * k3rx, y[i]);
		k4rx = vx[i] + dt * k3vx;

		vx[i+1] = vx[i] + (dt/6) * (k1vx + 2*k2vx + 2*k3vx + k4vx);
		x[i+1] = x[i] + (dt/6) * (k1rx + 2*k2rx + 2*k3rx + k4rx);


		k1vy = fy(x[i],y[i]);
		k1ry = vy[i];
		k2vy = fy(x[i], y[i] + (dt/2) * k1ry);
		k2ry = vy[i] + (dt/2) * k1vy;
		k3vy = fy(x[i], y[i] + (dt/2) * k2ry);
		k3ry = vy[i] + (dt/2) * k2vy;
		k4vy = fy(x[i], y[i] + dt * k3ry);
		k4ry = vy[i] + dt * k3vy;

		vy[i+1] = vy[i] + (dt/6) * (k1vy + 2*k2vy + 2*k3vy + k4vy);
		y[i+1] = y[i] + (dt/6) * (k1ry + 2*k2ry + 2*k3ry + k4ry);
	}

	// FLOp:
	//
	// (10 + 13 + 3 + 13 + 3 + 12 + 2 + 16) * 2 * N
	// 144 N
}


void semiI (double * x, double * y, double * vx, double * vy, double dt){

	for(size_t i = 0; i < NSTEPS; i++){
	vx[i+1] = vx[i] + fx(x[i],y[i]) * dt;
        vy[i+1] = vy[i] + fy(x[i],y[i]) * dt;

	x[i+1] = x[i] + vx[i+1] * dt;
        y[i+1] = y[i] + vy[i+1] * dt;
	}

	// FLOp:
	//
	// (12 * 2 + 2 * 2) N
	// 28 N


}

void leapfrog (double * x, double * y, double * vx, double * vy, double dt){

	for(size_t i = 0; i < NSTEPS; i++){
		x[i+1] = x[i] + vx[i] * dt + vx[i] * 0.5 * dt * dt;
		y[i+1] = y[i] + vy[i] * dt + vy[i] * 0.5 * dt * dt;
		vx[i+1] = vx[i] + (fx(x[i],y[i]) + fx(x[i+1],y[i+1]))*0.5*dt;
		vy[i+1] = vy[i] + (fy(x[i],y[i]) + fy(x[i+1],y[i+1]))*0.5*dt;
	}

	// FLOp:
	//
	// 6 * 2 + 24 * 2
	// 60 N

}

void init (double* &x, double* &y, double* &vx, double* &vy){

	// Add initial conditions
	x[0] = 1.0;
	y[0] = 0.0;
	vx[0] = 0.0;
	vy[0] = std::sqrt(1+E);


}

void write(double* x, double* y, double* vx, double* vy, size_t nsteps ,char* filename){


  std::ofstream outfile ((std::string) folder+"/" + (std::string) (filename));
  std::cout << "Path: " << (std::string) folder +"/" + (std::string) (filename) << "\n";

  if (outfile.is_open())
  {
    outfile << "nsteps X Y Vx Vy\n";

    for(size_t i = 0; i < nsteps; i++){
        outfile << nsteps <<" " << x[i] << " " << y[i] << " " << vx[i] << " " << vy[i] << "\n";
    }
    outfile.close();
  }
  else std::cout << "Unable to open file";


}

void check(double* x, double* y, double* vx, double* vy )
{

	for(size_t i = 0; i < NSTEPS; i++){

	if (std::sqrt(1 + E) - std::abs(x[i]*vy[i] - y[i]*vx[i]) > 1e-10) {

		std::cout << "Angular momentum conservation failed!\n";
		std::cout << std::sqrt(1 + E) - std::abs(x[i]*vx[i] - y[i]*vy[i] )<< "\n";
	}

	if (std::abs(0.5 * (1+E) - 1 / (std::sqrt(x[i]*x[i] + y[i]*y[i]))) - std::abs(-(1+E)/(1-E)) > 1e-10) {

		std::cout << "Energy conservation failed!\n";
	}
	}

}

int main(){


	const double dt =  T()/NSTEPS;

	double *x = (double*) malloc(NSTEPS*sizeof(double));
	double *y = (double*) malloc(NSTEPS*sizeof(double));
	double *vx = (double*) malloc(NSTEPS*sizeof(double));
	double *vy = (double*) malloc(NSTEPS*sizeof(double));
	double *times = (double*) malloc(5*sizeof(double));

	for(int i = 0; i < NSTEPS; i++){
	x[i] = 0;
	y[i] = 0;
	vx[i] = 0;
	vy[i] = 0;
	}

	void (*p[5]) (double* x, double* y, double* vx, double* vy, double dt);
	char* funNames[5] = {(char*) "eEuler.txt",(char*) "rk2.txt",(char*) "rk4.txt",(char*) "semiI.txt",(char*) "leapfrog.txt"};


	p[0] = eEuler;
	p[1] = rk2;
	p[2] = rk4;
	p[3] = semiI;
	p[4] = leapfrog;




	for(size_t i = 0; i < 5; i++){
		init(x,y,vx,vy);
		auto tStart = std::chrono::high_resolution_clock::now();
		(*p[i]) (x,y,vx,vy,dt);
		auto tEnd = std::chrono::high_resolution_clock::now();
		check(x,y,vx,vy);
	        times[i] = std::chrono::duration<double>(tEnd - tStart).count();
		write(x, y, vx, vy, NSTEPS, (char*)funNames[i]);

	}

	std::cout << "----------------------\n";
	std::cout << "Exec times\n";
	std::cout <<"E: " << times[0] << "\n";
	std::cout <<"RK2: " << times[1] << "\n";
	std::cout <<"RK4: " << times[2] << "\n";
	std::cout <<"SE: " << times[3] << "\n";
	std::cout <<"L: " << times[4] << "\n";
	std::cout << "----------------------\n";
	std::cout << "GFLOps/s\n";
	std::cout <<"E: " << 28 * NSTEPS / times[0] * 1e-9 << "\n";
	std::cout <<"RK2: " << 54 * NSTEPS / times[1] * 1e-9 << "\n";
	std::cout <<"RK4: " << 144 * NSTEPS / times[2] * 1e-9 << "\n";
	std::cout <<"SE: " << 28 * NSTEPS / times[3] * 1e-9 << "\n";
	std::cout <<"L: " << 60 * NSTEPS / times[4] * 1e-9 << "\n";
	std::cout << "----------------------\n";

	free(x);
	free(y);
	free(vx);
	free(vy);

return 0;
}
