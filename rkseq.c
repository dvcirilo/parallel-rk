/* Runge Kutta for a set of first order differential equations */

#include <stdio.h>
#include <math.h>

#define N 1                     /* number of first order equations */
#define dist 0.2                /* stepsize in t*/
#define MAX 5.0                /* max for t */

FILE *output;                   /* internal filename */

void runge4(double x, double y[], double step); /* Runge-Kutta function */

double F(double x, double y[], int i);          /* function for derivatives */

void runge4(double x, double y[], double step){
  double h=step/2.0,                      /* the midpoint */
          g1[N], g2[N], g3[N],            /* temporary storage arrays */
          g4[N], k1[N];      /* for Runge-Kutta */
  int i;

  for (i=0; i<N; i++){
    k1[i] = h*F(x+0.5*h, y, i);
    g1[i] = y[i] + 0.5 * k1[i];
    g2[i] = y[i] + (0.33) * k1[i];
    g3[i] = g2[i] + (0.33)*h*F(x+0.5*h, g2, i);
    g4[i] = F(x+(0.33)*h, g2, i);
    y[i] += h*(-2*F(x+0.5*h,g1,i)+(1.5)*g4[i]+(1.5)*F(x+(0.66)*h,g3,i));
  }
}

double F(double x, double y[], int i){
  if (i==0) return(1 - x*y[0]);                 /* derivative of first equation */
  //if (i==1) return(-0.2*y[1]-y[0]);       /* derivative of second equation */
}

int main(){
  double t, y[N];
  int j;

  output=fopen("data/rkseq.dat", "w");                   /* external filename */

  y[0]=1.0;                                       /* initial position */
  //y[1]=0.0;                                       /* initial velocity */
  fprintf(output, "0\t%f\n", y[0]);

  for (j=1; j*dist<=MAX ;j++)                     /* time loop */
  {
     t=j*dist;
     runge4(t, y, dist);
     fprintf(output, "%f\t%f\n", t, y[0]);
  }

  fclose(output);
}
