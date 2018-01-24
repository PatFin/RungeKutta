#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gnuplot_i.h" /* credits to N. Devillard */

/*------------------------ INITIAL VALUES AND CONSTANTS ----------------------*/
int STEPS_NUMBER = 2500;
double TIME_STEP = 0.01;

double INITIAL_U1_VALUE = 0.0;
double INITIAL_U2_VALUE = 0.0;
double INITIAL_T_VALUE = 0.0;

double A = 1.0; /* Parameter in the differential equation */
/*---------------------------- SHARED VARIABLES ------------------------------*/
double ** U;
double * T;
gnuplot_ctrl * gnuplot;


/*------------------------------- FUNCTIONS ----------------------------------*/
/* Fonction f for the computation of u1 */
double f1 (double t, double u1, double u2)
{
	/* In our particular situation, f1 only returns u2. Parameters t and u2 are 
	 * not used, this is not a mistake
	 */
	return u2;
}

/* Fonction f for the computation of u2 */
double f2 (double t, double u1, double u2)
{
	/* Parameter t is not used in our particular case, this is not an error */
	return 1 - A*u2 - u1;
}

/* Routine that sets up the initial conditions of the system */
void init ()
{
  U = malloc( sizeof(double *) * 2 );
	U [0] = malloc( sizeof(double) * STEPS_NUMBER );
	U [1] = malloc( sizeof(double) * STEPS_NUMBER );
	T = malloc( sizeof(double) * STEPS_NUMBER );

  U [0][ 0 ] = INITIAL_U1_VALUE;
	U [1][ 0 ] = INITIAL_U2_VALUE;
	T [ 0 ] = INITIAL_T_VALUE;
  gnuplot = gnuplot_init();
}

/* Computes the "lastStep + 1" approximation and stores it into the X array */
void nextStep (int lastStepIndex)
{
  /* The next approximation of U (U[1] and U[2] is computed using the 
	 * Runge-Kutta approach.
	 */
  double k1, k2, k3, k4, t, u1, u2;
	/* Retrieving the variables at the nth approximation */
	u1 = U [0][ lastStepIndex ];
	u2 = U [1][ lastStepIndex ];
	t = T [ lastStepIndex ];

	/* Computing the next u1 */
  k1 = f1 (t, u1, u2);
  k2 = f1 (t + (TIME_STEP/2) , u1 + (TIME_STEP*k1/2), u2 + (TIME_STEP*k1/2) );
  k3 = f1 (t + (TIME_STEP/2), u1 + (TIME_STEP*k2/2), u2 + (TIME_STEP*k2/2) );
  k4 = f1 (t + TIME_STEP, u1 + (TIME_STEP*k3), u2 + (TIME_STEP*k3) );

  U [0][ lastStepIndex + 1 ] = u1 + (TIME_STEP/6)*(k1 + k2 + k3 + k4);
	
	/* Computing the next u2, we recycle variables k1 to k4 */
	k1 = f2 (t, u1, u2);
  k2 = f2 (t + (TIME_STEP/2) , u1 + (TIME_STEP*k1/2), u2 + (TIME_STEP*k1/2) );
  k3 = f2 (t + (TIME_STEP/2), u1 + (TIME_STEP*k2/2), u2 + (TIME_STEP*k2/2) );
  k4 = f2 (t + TIME_STEP, u1 + (TIME_STEP*k3), u2 + (TIME_STEP*k3) );

  U [1][ lastStepIndex + 1 ] = u2 + (TIME_STEP/6)*(k1 + k2 + k3 + k4);

	T [ lastStepIndex + 1 ] = t + TIME_STEP;
}

/* Clears allocated memory and closes opened files */
void end ()
{
  free(U [0]);
	free(U [1]);
	free(U);
	free(T);
	gnuplot_close(gnuplot);
}

/*---------------------------------- MAIN ------------------------------------*/

/* The user can adjust the value of the parameter a of the studied PDE. 
 * By default it is one but if something is written after the name of the
 * executable it will be changed to 10.
 */
int main (int argc, char * args)
{
	/* Adapting the parameters for the equation and the display */
	if ( argc > 1 )
	{
		A = 10;
		STEPS_NUMBER = 8000;
	}


	/* Computation */
  init();

  int i;
  for (i = 0; i < STEPS_NUMBER - 1 ; i++)
  {
    nextStep(i);
    printf("%f %f %f\n", T[i], U[0][i], U[1][i]);
  }

	/* Preparing the display */
	char * legend = "A = 1";
	if ( argc > 1 )
	{
		legend = "A = 10";
	}
	/* Display */
	gnuplot_plot_xy(gnuplot, T, U [1], STEPS_NUMBER, legend);
	sleep(10);

	/* Cleanup and exit */
	end();
  return 0;
}
