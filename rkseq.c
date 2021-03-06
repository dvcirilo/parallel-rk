/* Runge Kutta for a set of first order differential equations */

#include <stdio.h>
#include <math.h>

#define dist 0.2                /* stepsize in t*/
#define MAX 5.0                /* max for t */

FILE *output;                   /* internal filename */

double runge4(double t, double y); /* Runge-Kutta function */

double F(double x, double y);          /* function for derivatives */

double runge4(double t, double y){
  double h = dist, k1, k2, k3, k4;
  int i;

  k1 = h*F(t, y);
  k2 = h*F(t + 0.5*h, y+0.5*k1);
  k3 = h*F(t + 0.5*h, y+0.5*k2);
  k4 = h*F(t + h, y + k3);

  y += 0.33*(k1 + 2*k2 + 2*k3 + k4);
  return y;
}

double F(double x, double y){
  return(1 - x*y);
}

int main(){
  double t, y;
  int j;

  output=fopen("data/rkseq.dat", "w");                   /* external filename */

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
