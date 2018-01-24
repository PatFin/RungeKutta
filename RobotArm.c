#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "gnuplot_i.h" /* credits to N. Devillard */

/*------------------------ INITIAL VALUES AND CONSTANTS ----------------------*/
int STEPS_NUMBER = 150;
double TIME_STEP = 0.1;

const double INITIAL_V_VALUE = 0.0;
const double INITIAL_W_VALUE = 0.0;
const double INITIAL_X_VALUE = 0.0;
const double INITIAL_Y_VALUE = 0.0;
const double INITIAL_T_VALUE = 0.0;

const double TARGET_V = 1.6;	/* Target V angle in radians */
const double TARGET_W = 0.7;  /* Target W angle in radians */

const double GRAVITY = 9.81;	/* Earth acceleration */

const double L1 = 1.0;	/* Length of the first rod */
const double R1 = 0.5;	/* Distance from its anchor point to center of mass */
const double M1 = 1.0;	/* Mass of the first rod */
const double IZ1 = 5.0;	/* Moment of inertia of the first rod */

const double L2 = 1.0;	/* Length of the second rod */
const double R2 = 0.5;	/* Distance from its anchor point to center of mass */
const double M2 = 1.0;	/* Mass of the second rod */
const double IZ2 = 5.0;	/* Moment of inertia of the second rod */

const double Kp1 = 5.0;	/* Parameter of the input wrt current angle for v */
const double Ku1 = 15.0;/* Parameter of the input wrt angular velocity for v */
const double Kp2 = 10.0;/* Parameter of the input wrt current angle for w */
const double Ku2 = 10.0;/* Parameter of the input wrt angular velocity for w */

double ALPHA, BETA, DELTA;

double * Xcords;	/* Set of X coordinates to display */
double * Ycords;	/* Set of Y coordinates to display */
const int MAX_DISPLAY = 100;
/*---------------------------- SHARED VARIABLES ------------------------------*/
double V, W, X, Y, T;
gnuplot_ctrl * gnuplot;


/*------------------------------- FUNCTIONS ----------------------------------*/
double t1 (double t, double v, double w, double x, double y)
{
	return Kp1*(TARGET_V - v) - Ku1*x;
}

double t2 (double t, double v, double w, double x, double y)
{
	return Kp2*(TARGET_W - w) - Ku2*y;
}

double c1 (double t, double v, double w, double x, double y)
{
	return BETA*sin(w)*((x + y)*y - y*x);
}

double c2 (double t, double v, double w, double x, double y)
{
	return BETA*sin(w)*x*x;
}

double g1 (double t, double v, double w, double x, double y)
{
	return 0;
	return -GRAVITY*(R1*M1*cos(v) + (L1*cos(v) + R2*cos(w))*M2);
}

double g2 (double t, double v, double w, double x, double y)
{
	return 0;
	return -GRAVITY*(R2*M2*cos(w));
}

double m1 (double t, double v, double w, double x, double y)
{
	return (ALPHA + 2*BETA*cos(w)) + (DELTA + BETA*cos(w));
}

double m2 (double t, double v, double w, double x, double y)
{
	return (DELTA + BETA*cos(w)) + DELTA;
}


/* Fonction f for the computation of v */
double fv (double t, double v, double w, double x, double y)
{
	return x;
}

/* Fonction f for the computation of w */
double fw (double t, double v, double w, double x, double y)
{
	return y;
}

/* Fonction f for the computation of x */
double fx (double t, double v, double w, double x, double y)
{
	return (t1(t,v,w,x,y) - c1(t,v,w,x,y) - g1(t,v,w,x,y))/m1(t,v,w,x,y);
}

/* Fonction f for the computation of y */
double fy (double t, double v, double w, double x, double y)
{
	return (t2(t,v,w,x,y) - c2(t,v,w,x,y) - g2(t,v,w,x,y))/m2(t,v,w,x,y);
}

/* Routine that sets up the initial conditions of the system */
void init ()
{
	ALPHA = IZ1 + IZ2 + M1*R1*R1 + M2*(L1*L1 + R2*R2);
	BETA = M2*L1*R2;
	DELTA = IZ2 + M2*R2*R2;

	V = INITIAL_V_VALUE;
	W = INITIAL_W_VALUE;
	X = INITIAL_X_VALUE;
	Y = INITIAL_Y_VALUE;
	T = INITIAL_T_VALUE;
	gnuplot = gnuplot_init();
	gnuplot_cmd(gnuplot, "set yr [-2:2]");
	gnuplot_cmd(gnuplot, "set xr [-2:2]");
	
	Xcords = malloc( sizeof(double) * MAX_DISPLAY);
	Ycords = malloc( sizeof(double) * MAX_DISPLAY);
}

/* Computes the next approximation of the system and updates the values in the 
 * variables V, W, X, Y and T.
 */
void nextStep ()
{
	/* The next approximations of [v,w,x,y] are computed using the 
	 * Runge-Kutta approach.
	 */
	double k1, k2, k3, k4, t, v, w, x, y;
	/* Retrieving the variables at the nth approximation */
	v = V;
	w = W;
	x = X;
	y = Y;
	t = T;

	/* Computing the next v */
	k1 = fv (t, v, w, x ,y);
	k2 = fv (t + (TIME_STEP/2) , v + (TIME_STEP*k1/2), w + (TIME_STEP*k1/2), 
			x + (TIME_STEP*k1/2), y + (TIME_STEP*k1/2) );
	k3 = fv (t + (TIME_STEP/2), v + (TIME_STEP*k2/2), w + (TIME_STEP*k2/2), 
			x + (TIME_STEP*k2/2), y + (TIME_STEP*k2/2) );
	k4 = fv (t + TIME_STEP, v + (TIME_STEP*k3), w + (TIME_STEP*k3),
			x + (TIME_STEP*k3), y + (TIME_STEP*k3) );

	V = v + (TIME_STEP/6)*(k1 + k2 + k3 + k4);
	
	/* Computing the next w, we recycle variables k1 to k4 */
	k1 = fw (t, v, w, x ,y);
	k2 = fw (t + (TIME_STEP/2) , v + (TIME_STEP*k1/2), w + (TIME_STEP*k1/2), 
			x + (TIME_STEP*k1/2), y + (TIME_STEP*k1/2) );
	k3 = fw (t + (TIME_STEP/2), v + (TIME_STEP*k2/2), w + (TIME_STEP*k2/2), 
			x + (TIME_STEP*k2/2), y + (TIME_STEP*k2/2) );
	k4 = fw (t + TIME_STEP, v + (TIME_STEP*k3), w + (TIME_STEP*k3),
			x + (TIME_STEP*k3), y + (TIME_STEP*k3) );

	W = w + (TIME_STEP/6)*(k1 + k2 + k3 + k4);

	/* Computing the next x, we recycle variables k1 to k4 */
	k1 = fx (t, v, w, x ,y);
	k2 = fx (t + (TIME_STEP/2) , v + (TIME_STEP*k1/2), w + (TIME_STEP*k1/2), 
			x + (TIME_STEP*k1/2), y + (TIME_STEP*k1/2) );
	k3 = fx (t + (TIME_STEP/2), v + (TIME_STEP*k2/2), w + (TIME_STEP*k2/2), 
			x + (TIME_STEP*k2/2), y + (TIME_STEP*k2/2) );
	k4 = fx (t + TIME_STEP, v + (TIME_STEP*k3), w + (TIME_STEP*k3),
			x + (TIME_STEP*k3), y + (TIME_STEP*k3) );

	X = x + (TIME_STEP/6)*(k1 + k2 + k3 + k4);

	/* Computing the next w, we recycle variables k1 to k4 */
	k1 = fy (t, v, w, x ,y);
	k2 = fy (t + (TIME_STEP/2) , v + (TIME_STEP*k1/2), w + (TIME_STEP*k1/2), 
			x + (TIME_STEP*k1/2), y + (TIME_STEP*k1/2) );
	k3 = fy (t + (TIME_STEP/2), v + (TIME_STEP*k2/2), w + (TIME_STEP*k2/2), 
			x + (TIME_STEP*k2/2), y + (TIME_STEP*k2/2) );
	k4 = fy (t + TIME_STEP, v + (TIME_STEP*k3), w + (TIME_STEP*k3),
			x + (TIME_STEP*k3), y + (TIME_STEP*k3) );

	Y = y + (TIME_STEP/6)*(k1 + k2 + k3 + k4);

	T = t + TIME_STEP;
}

/* Computes the new positions of the robot arms and displays it */
void displayPoints ()
{
	static int dispIndex = 0;
	static int toDisplay = 0;
	
	/* Retrieving the current angles of the system */
	double theta1 = V;
	double theta2 = W;
	
	/* Computing the coordinates of the first point */
	double x1, y1;
	x1 = cos(theta1) * L1;
	y1 = sin(theta1) * L1;
	
	/* Computing the coordinates of the second point */
	double x2, y2;
	x2 = x1 + cos(theta2)*L2;
	y2 = y1 + sin(theta2)*L2;
	
	Xcords [dispIndex] = x1;
	Ycords [dispIndex] = y1;
	
	dispIndex = (dispIndex+1)%MAX_DISPLAY;
	
	Xcords [dispIndex] = x2;
	Ycords [dispIndex] = y2;
	
	dispIndex = (dispIndex+1)%MAX_DISPLAY;
	
	toDisplay = toDisplay + 2;
	if (toDisplay > MAX_DISPLAY)
	{
		toDisplay = MAX_DISPLAY;
	}
	gnuplot_resetplot(gnuplot);
	gnuplot_plot_xy(gnuplot, Xcords, Ycords, toDisplay, ""); 
}

/* Clears allocated memory and closes opened files */
void end ()
{
	free(Xcords);
	free(Ycords);
	gnuplot_close(gnuplot);
}

/*---------------------------------- MAIN ------------------------------------*/

/* The user can adjust the value of the parameter a of the studied PDE. 
 * By default it is one but if something is written after the name of the
 * executable it will be changed to 10.
 */
int main (int argc, char * args)
{
	/* Computation */
	init();

	printf("%f %f %f %f %f\n", T, V, W, X, Y);

	int i;
	for (i = 0; i < STEPS_NUMBER - 1 ; i++)
	{
		nextStep();
		printf("%i %f %f %f %f %f\n", i, T, V, W, X, Y);
		
		displayPoints();
		sleep(1);
	}

	
	/* Display 
	gnuplot_plot_xy(gnuplot, T, U [1], STEPS_NUMBER, "ROBOT ARM");
	sleep(10); */

	/* Cleanup and exit */
	end();
  return 0;
}
