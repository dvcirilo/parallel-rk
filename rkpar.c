#include <stdio.h>
#include <math.h>

#define dist 0.2                /* stepsize in t*/
#define MAX 5.0                /* max for t */

FILE *output;                   /* internal filename */

double runge4(double t, double y); /* Runge-Kutta function */

double F(double x, double y);          /* function for derivatives */

double runge4(double t, double y){
  double h = dist, hh = h/2, h3 = h/3, k1, k2, k3, k4;
  int i;

  k1 = y;
  k2 = y + hh*F(t + hh, k1);
  k3 = y + h3*F(t + hh, k1);
  k4 = y + h3*F(t + hh, k1) + h3*F(t + hh, k3);

  y += h*(-2*F(t + hh, k2) + 1.5*F(t + h3, k3) + 1.5*F(t + 2*h3, k4));
  return y;
}

double F(double x, double y){
  return(1 - x*y);
}

int main(){
  double t, y;
  int j;

  output=fopen("data/rkpar.dat", "w");                   /* external filename */

  y = 1.0;                                       /* initial position */
  fprintf(output, "0\t%f\n", y);

  for (j=1; j*dist<=MAX ;j++)                     /* time loop */
  {
     t=j*dist;
     y = runge4(t, y);
     fprintf(output, "%f\t%f\n", t, y);
  }

  fclose(output);
}
