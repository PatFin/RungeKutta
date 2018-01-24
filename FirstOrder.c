#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gnuplot_i.h" /* credits to N. Devillard */

/*------------------------ INITIAL VALUES AND CONSTANTS ----------------------*/
int STEPS_NUMBER = 1000;
double TIME_STEP = 0.01;

double INITIAL_X_VALUE = 0.0;
double INITIAL_T_VALUE = 0.0;

/*---------------------------- SHARED VARIABLES ------------------------------*/
double * X;
double * T;
gnuplot_ctrl * gnuplot;


/*------------------------------- FUNCTIONS ----------------------------------*/
/* Fonction f as described in the Runge-Kutta method */
double f (double t, double x)
{
	/* In our particular situation, time 't' does not intervene in the 
	 * formula. This is not a mistake 
	*/
	return -x+1;
}

/* Routine that sets up the initial conditions of the system */
void init ()
{
  X = malloc( sizeof(double) * STEPS_NUMBER );
	T = malloc( sizeof(double) * STEPS_NUMBER );

  X [ 0 ] = INITIAL_X_VALUE;
	T [ 0 ] = INITIAL_T_VALUE;
  gnuplot = gnuplot_init();
}

/* Computes the "lastStep + 1" approximation and stores it into the X array */
void nextStep (int lastStepIndex)
{
  /* The next approximation of X is computed using the Runge-Kutta approach */
  double k1, k2, k3, k4, t, x;
	x = X [ lastStepIndex ];
	t = T [ lastStepIndex ];
  k1 = f (t, x);
  k2 = f (t + (TIME_STEP/2) , x + (TIME_STEP*k1/2) );
  k3 = f (t + (TIME_STEP/2), x + (TIME_STEP*k2/2) );
  k4 = f (t + TIME_STEP, x + (TIME_STEP*k3) );

  X [ lastStepIndex + 1 ] = x + (TIME_STEP/6)*(k1 + k2 + k3 + k4);
	T [ lastStepIndex + 1 ] = t + TIME_STEP;
}

/* Clears allocated memory and closes opened files */
void end ()
{
  free(X);
	free(T);
	gnuplot_close(gnuplot);
}

/*----------------------------------- MAIN -----------------------------------*/
int main ()
{
  init();

  int i;
  for ( i = 0; i < STEPS_NUMBER - 1 ; i ++)
  {
    nextStep(i);
    printf("%f %f\n", T[i], X[i] );
  }
	gnuplot_plot_xy(gnuplot, T, X, STEPS_NUMBER, "Exercise 1");

	sleep(10);


	end();
  return 0;
}
