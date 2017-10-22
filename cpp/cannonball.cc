#include <cmath>
#include <iostream>

using namespace std;

// Constants
const double m = 2.0;
const double g = 9.81;
const double dt = 0.1;
const double total_t = 100;
const int steps = total_t/dt;
  


void computeForces( double f[]){
  f[0] = 0.0;
  f[1] = -m*g;
}

void stepEuler( double x[], double v[], double f[]){
  computeForces(f);
  x[0] += v[0]*dt;
  x[1] += v[1]*dt;
  v[0] += f[0]*dt/m;
  v[1] += f[1]*dt/m;
}

double energy( double x[], double v[]){
  double energy = 0.5*m*(v[0]*v[0]+v[1]*v[1])+m*g*x[1];
  return energy;
}

int main(){

  // Initialize Variables
  double t = 0.0;
  double x[2] = {0.0, 0.0};
  double v[2] = {50.0, 50.0};
  double f[2] = {0.0, 0.0};
  double traj[steps][2];
  
  while(t <= total_t){
    cout << x[0] << " " << x[1] << " " << energy(x,v) << endl;
    stepEuler(x, v, f);
    t = t+dt;
  }
  return 0;
}
