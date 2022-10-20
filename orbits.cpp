#include<iostream>
#include<cmath>
#include<omp.h>
#include<fstream>
#include<chrono>

# define E 0.9
# define NSTEPS 1e5

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

		vx[i+1] = vx[i] + fx(x[i+1],y[i+1]) * dt;
		vy[i+1] = vy[i] + fy(x[i+1],y[i+1]) * dt;
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


  std::ofstream outfile (filename);
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

int main(){
	const double dt =  T()/NSTEPS;

	double *x = (double*) malloc(NSTEPS*sizeof(double));
	double *y = (double*) malloc(NSTEPS*sizeof(double));
	double *vx = (double*) malloc(NSTEPS*sizeof(double));
	double *vy = (double*) malloc(NSTEPS*sizeof(double));

	for(size_t i = 0; i < NSTEPS; i++){
	x[i] = 0;
	y[i] = 0;
	vx[i] = 0;
	vy[i] = 0;
	}

	// Explicit Euler

	init(x,y,vx,vy);
	auto tStartE = std::chrono::high_resolution_clock::now();
	eEuler(x,y,vx,vy,dt);
	auto tEndE = std::chrono::high_resolution_clock::now();
	write(x,y,vx,vy,NSTEPS,(char*) "eEuler.txt");

	auto tExecE = std::chrono::duration<double>(tEndE - tStartE).count();

	// Runge Kutta 2
	init(x,y,vx,vy);
	auto tStartRK2 = std::chrono::high_resolution_clock::now();
	rk2(x,y,vx,vy,dt);
	auto tEndRK2 = std::chrono::high_resolution_clock::now();

	write(x,y,vx,vy,NSTEPS,(char*) "rk2.txt");
	auto tExecRK2 = std::chrono::duration<double>(tEndRK2 - tStartRK2).count();

	// Runge Kutta 4
	init(x,y,vx,vy);
	auto tStartRK4 = std::chrono::high_resolution_clock::now();
	rk4(x,y,vx,vy,dt);
	auto tEndRK4 = std::chrono::high_resolution_clock::now();
	write(x,y,vx,vy,NSTEPS,(char*) "rk4.txt");
	auto tExecRK4 = std::chrono::duration<double>(tEndRK4 - tStartRK4).count();

	// Semi Implicit Euler
	init(x,y,vx,vy);
	auto tStartSE = std::chrono::high_resolution_clock::now();
	semiI(x,y,vx,vy,dt);
	auto tEndSE = std::chrono::high_resolution_clock::now();
	write(x,y,vx,vy,NSTEPS,(char*) "semiI.txt");
	auto tExecSE = std::chrono::duration<double>(tEndSE - tStartSE).count();

	// Leapfrog
	init(x,y,vx,vy);
	auto tStartL = std::chrono::high_resolution_clock::now();
	leapfrog(x,y,vx,vy,dt);
	auto tEndL = std::chrono::high_resolution_clock::now();
	write(x,y,vx,vy,NSTEPS,(char*) "leapfrog.txt");
	auto tExecL = std::chrono::duration<double>(tEndL - tStartL).count();
	//28,54,144,28,60
	/*
	std::cout << "Exec times\n";
	std::cout <<"E: " <<tExecE << "\n";
	std::cout <<"RK2: " << tExecRK2 << "\n";
	std::cout <<"RK4: " << tExecRK4 << "\n";
	std::cout <<"SE: " << tExecSE << "\n";
	std::cout <<"L: " <<tExecL << "\n";
	*/
	std::cout << "GFLOps/s\n";
	std::cout <<"E: " << 28 * NSTEPS / tExecE * 1e-9 << "\n";
	std::cout <<"RK2: " << 54 * NSTEPS / tExecRK2 * 1e-9 << "\n";
	std::cout <<"RK4: " << 144 * NSTEPS / tExecRK4 * 1e-9 << "\n";
	std::cout <<"SE: " << 28 * NSTEPS / tExecSE * 1e-9 << "\n";
	std::cout <<"L: " << 60 * NSTEPS / tExecL * 1e-9 << "\n";

	free(x);
	free(y);
	free(vx);
	free(vy);


return 0;
}
