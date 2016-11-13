
////////////////////////////////////////////////////////////////////////////////
// File: embedded_verner_8_9.c                                                //
// Routines:                                                                  //
//    Embedded_Verner_8_9                                                     //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     The Runge-Kutta-Verner method is an adaptive procedure for approxi-    //
//     mating the solution of the differential equation y'(x) = f(x,y) with   //
//     initial condition y(x0) = c.  This implementation evaluates f(x,y)     //
//     sixteen times per step using embedded eighth order and ninth order     //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * ( 103/1680 * k1 - 27/140 * k8 + 76/105 * k9    //
//              - 201/280 * k10 + 1024/1365 * k11 + 3/7280 k12 + 12/35 k13    //
//                                                             + 9/280 k14 )  //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/12, y[i] + h*k1/12),                                    //
//     k3 = f( x[i]+h/9, y[i]+h/27*( k1 + 2 k2) ),                            //
//     k4 = f( x[i]+h/6, y[i]+h/24*( k1 + 3 k3) ),                            //
//     k5 = f( x[i]+(2+2s6)h/15, y[i]+h/375*((4+94s6) k1 - (282+252s6) k3     //
//                                                      + (328+208s6) k4)),   //
//     k6 = f( x[i]+(6+s6)h/15, y[i]+h*( (9-s6) k1/150 + (312+32s6) k4/1425   //
//                                                   + (69+29s6) k5/570 ) ),  //
//     k7 = f( x[i]+(6-s6)h/15, y[i]+h*( (927-347s6) k1/ 1250 +               //
//                          (-16248+7328s6) k4/9375 + (-489+179s6) k5 /3750   //
//                                               + (14268-5798s6) k6/9375) ), //
//     k8 = f( x[i]+2h/3, y[i]+h/54*( 4 k1 + (16-s6) k6 + (16+s6) k7) )       //
//     k9 = f( x[i]+h/2, y[i]+h/512*( 38 k1 + (118-23s6) k6 + (118+23s6) k7   //
//                                                                - 18 k8) ), //
//     k10 = f( x[i]+h/3, y[i]+h*( 11 k1/144 + (266-s6) k6/864                //
//                                   + (266+s6) k7/864 - K8/16 - 8 k9/27 ) )  //
//     k11 = f( x[i]+h/4, y[i]+h*( (5034-271s6) k1/61440                      //
//                       + (7859-1626s6) k7/10240 + (-2232+813s6) k8/20480 +  //
//                             (-594+271s6) k9/960 + (657-813s6) k10/5120) )  //
//     k12 = f( x[i]+4h/3, y[i]+h*( (5996-3794s6) k1/405 + (-4342-338s6) k6/9 //
//                        + (154922-40458s6)k7/135 + (-4176+3794s6) k8/45     //
//                        + (-340864+242816s6)k9/405 + (26304-15176)s6 k10/45 //
//                                                        - 26624 k11 /81     //
//     k13 = f( x[i]+5h/6, y[i]+h*( (3793 + 2168s6) k1 / 103680               //
//                 + (4042+2263s6) k6 / 13824 + (-231278+40717s6) k7 / 69120  //
//                      + (7947 - 2168s6) k8 / 11520 + (1048-542s6) k9 / 405  //
//              + (-1383+542s6) k10 / 720 + 2624 k11 / 1053 + 3 k12 / 1664 )  //
//     k14 = f( x[i]+h, y[i]+h*( - 137 k1 / 1296 + (5642-337s6) k6 / 864      //
//         + (5642+337s6) k7 / 864 - 299 k8 / 48 + 184 k9 / 81 - 44 k10 / 9   //
//                          - 5120 k11 / 1053 - 11 k12 / 468 + 16 k13 / 9 ),  //
//     k15 = f( x[i]+h/6, y[i]+h*( (33617 - 2168s6) k1 / 518400               //
//                 + (-3846+31s6) k6 / 13824 + (155338-52807s6) k7 / 345600   //
//                      + (-12537 + 2168s6 k8 / 57600 + (92+542s6) k9 / 2025  //
//               + (-1797-542s6) k10 / 3600 + 320 k11 / 567 - k12 / 1920      //
//                                                        + 4 k13 / 105 ),    //
//     k16 = f( x[i]+h, y[i]+h*( (-36487 - 30352s6) k1 / 279600               //
//            + (-29666-4499s6) k6 / 7456 + (2779182-615973s6) k7 / 186400    //
//            + (-94329 + 91056s6 k8 / 93200 + (-232192+121408s6) k9 / 17475  //
//                      + (101226-22764s6) k10 / 5825 - 169984 k11 / 9087     //
//                      - 87 k12 / 30290 + 492 k13 / 1165 + 1260 k15 / 233 )  //
//     and s6 = sqrt(6), x[i+1] = x[i] + h.                                   //
//                                                                            //
//     The error is estimated to be                                           //
//        err = - h*( 1911 k1 - 34398 k8 + 61152 k9 - 114660 k10 + 114688 k11 //
//           + 63 k12 + 13104 k13 + 3510 k14 - 39312 k15 - 6058 k16) / 109200  //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/8     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
//                                                                            //
//  Reference:                                                                //
//     SIAM Journal on Numerical Analysis Vol 15 No. 4 Aug 1978 pages 772-790.//
//     "EXPLICIT RUNGE-KUTTA METHODS WITH ESTIMATES OF THE LOCAL TRUNCATION   //
//     ERROR" - J.H. VERNER.                                                  //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

#define ATTEMPTS 12
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0

static double Runge_Kutta(double (*f)(double,double), double *y, double x,
                                                                   double h);

////////////////////////////////////////////////////////////////////////////////
// int Embedded_Verner_8_9( double (*f)(double, double), double y[],          //
//       double x, double h, double xmax, double *h_next, double tolerance )  //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation y'=f(x,y) with the      //
//     initial condition y(x) = y[0].  The value at xmax is returned in y[1]. //
//     The function returns 0 if successful or -1 if it fails.                //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y0) corresponding to the //
//                initial condition y(x0) = y0.                               //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at xmax.                               //
//     double x   The initial value of x.                                     //
//     double h   Initial step size.                                          //
//     double xmax The endpoint of x.                                         //
//     double *h_next   A pointer to the estimated step size for successive   //
//                      calls to Embedded_Verner_8_9.                         //
//     double tolerance The tolerance of y(xmax), i.e. a solution is sought   //
//                so that the relative error < tolerance.                     //
//                                                                            //
//  Return Values:                                                            //
//     0   The solution of y' = f(x,y) from x to xmax is stored y[1] and      //
//         h_next has the value to the next size to try.                      //
//    -1   The solution of y' = f(x,y) from x to xmax failed.                 //
//    -2   Failed because either xmax < x or the step size h <= 0.            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Embedded_Verner_8_9( double (*f)(double, double), double y[], double x,
                   double h, double xmax, double *h_next, double tolerance ) {

   static const double err_exponent = 1.0 / 8.0;

   double scale;
   double temp_y[2];
   double err;
   double yy;
   int i;
   int last_interval = 0;

      // Verify that the step size is positive and that the upper endpoint //
      // of integration is greater than the initial endpoint.              //

   if (xmax < x || h <= 0.0) return -2;
   
       // If the upper endpoint of the independent variable agrees with the //
       // initial value of the independent variable.  Set the value of the  //
       // dependent variable and return success.                            //

   *h_next = h;
   y[1] = y[0];
   if (xmax == x) return 0; 

        // Redefine the error tolerance to an error tolerance per unit    //
        // length of the integration interval.                            //

   tolerance /= (xmax - x);

        // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
        // maintain an error less than tolerance * (xmax-x) using an     //
        // initial step size of h and initial value: y = y[0]            //

   temp_y[0] = y[0];
   while ( !last_interval ) {
      for (i = 0; i < ATTEMPTS; i++) {
         last_interval = 0;
         if ( x + h >= xmax ) { h = xmax - x; last_interval = 1; }
         else if ( x + h + 0.5 * h > xmax ) h = 0.5 * (xmax - x);
         err = fabs( Runge_Kutta(f, temp_y, x, h) );
         if (err == 0.0) { scale = MAX_SCALE_FACTOR; break; }
         yy = (temp_y[0] == 0.0) ? tolerance : fabs(temp_y[0]);
         scale = 0.8 * pow( tolerance * yy /  err , err_exponent );
         scale = min( max(scale,MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
         if ( err < ( tolerance * yy ) ) break;
         h *= scale;
      }
      if ( i >= ATTEMPTS ) { *h_next = h; return -1; };
      temp_y[0] = temp_y[1];         
      x += h;
   }
   *h_next = h;
   y[1] = temp_y[1];
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Verner's embedded 8th and 9th order methods to       //
//     approximate the solution of the differential equation y'=f(x,y) with   //
//     the initial condition y = y[0] at x = x0.  The value at x + h is       //
//     returned in y[1].  The function returns err / h ( the absolute error   //
//     per step size ).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y[0]).                   //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at x + h.                              //
//     double x   Initial value of x.                                         //
//     double h   Step size                                                   //
//                                                                            //
//  Return Values:                                                            //
//     This routine returns the err / h.  The solution of y(x) at x + h is    //
//     returned in y[1].                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static double Runge_Kutta(double (*f)(double,double), double *y, double x0,
                                                                   double h) {

   #define SQRT6 2.449489742783178098197284074705891

                  // c2 = c3 = ... = c7 = 0, c15 = c16 = 0 //

   static const double c1 = 103.0 / 1680.0;
   static const double c8 = -27.0 / 140.0;
   static const double c9 = 76.0 / 105.0;
   static const double c10 = -201.0 / 280.0;
   static const double c11 = 1024.0 / 1365.0;
   static const double c12 = 3.0 / 7280.0;
   static const double c13 = 12.0 / 35.0;
   static const double c14 = 9.0 / 280.0;

                          // a1 = 0, a14 = a16 = 1 //

   static const double a2 = 1.0 / 12.0;
   static const double a3 = 1.0 / 9.0;
   static const double a4 = 1.0 / 6.0;
   static const double a5 = 2.0 * (1.0 + SQRT6) / 15.0;
   static const double a6 = (6.0 + SQRT6) / 15.0;
   static const double a7 = (6.0 - SQRT6) / 15.0;
   static const double a8 = 2.0 / 3.0;
   static const double a9 = 1.0 / 2.0;
   static const double a10 = 1.0 / 3.0;
   static const double a11 = 1.0 / 4.0;
   static const double a12 = 4.0 / 3.0;
   static const double a13 = 5.0 / 6.0;
   static const double a15 = 1.0 / 6.0;

             // b21 = 1/12, remaining missing bij's j < i are 0 //

   static const double b31 = 1.0 / 27.0;
   static const double b32 = 2.0 / 27.0;
   static const double b41 = 1.0 / 24.0;
   static const double b43 = 3.0 / 24.0;
   static const double b51 = (4.0 + 94.0 * SQRT6) / 375.0;
   static const double b53 = -(282.0 + 252.0 * SQRT6) / 375.0;
   static const double b54 = (328.0 + 208.0 * SQRT6) / 375.0;
   static const double b61 = (9.0 - SQRT6) / 150.0;
   static const double b64 = (312.0 + 32.0 * SQRT6) / 1425.0;
   static const double b65 = (69.0 + 29.0 * SQRT6) / 570.0;
   static const double b71 = (927.0 - 347.0 * SQRT6) / 1250.0;
   static const double b74 = (-16248.0 + 7328.0 * SQRT6) / 9375.0;
   static const double b75 = (-489.0 + 179.0 * SQRT6) / 3750.0;
   static const double b76 = (14268.0 - 5798.0 * SQRT6) / 9375.0;
   static const double b81 = 4.0 / 54.0;
   static const double b86 = (16.0 - SQRT6) / 54.0;
   static const double b87 = (16.0 + SQRT6) / 54.0;
   static const double b91 = 38.0 / 512.0;
   static const double b96 = (118.0 - 23.0 * SQRT6) / 512.0;
   static const double b97 = (118.0 + 23.0 * SQRT6) / 512.0;
   static const double b98 = - 18.0 / 512.0;
   static const double b10_1 = 11.0 / 144.0;
   static const double b10_6 = (266.0 - SQRT6) / 864.0;
   static const double b10_7 = (266.0 + SQRT6) / 864.0;
   static const double b10_8 = - 1.0 / 16.0;
   static const double b10_9 = - 8.0 / 27.0;
   static const double b11_1 = (5034.0 - 271.0 * SQRT6) / 61440.0;
   static const double b11_7 = (7859.0 - 1626.0 * SQRT6) / 10240.0;
   static const double b11_8 = (-2232.0 + 813.0 * SQRT6) / 20480.0;
   static const double b11_9 = (-594.0  + 271.0 * SQRT6) / 960.0;
   static const double b11_10 = (657.0 - 813.0 * SQRT6) / 5120.0;
   static const double b12_1 = (5996.0 - 3794.0 * SQRT6) / 405.0;
   static const double b12_6 = (-4342.0 - 338.0 * SQRT6) / 9.0;
   static const double b12_7 = (154922.0 - 40458.0 * SQRT6) / 135.0;
   static const double b12_8 = (-4176.0 + 3794.0 * SQRT6) / 45.0;
   static const double b12_9 = (-340864.0 + 242816.0 * SQRT6) / 405.0;
   static const double b12_10 = (26304.0 - 15176.0 * SQRT6) / 45.0;
   static const double b12_11 = -26624.0 / 81.0;
   static const double b13_1 = (3793.0 + 2168.0 * SQRT6) / 103680.0;
   static const double b13_6 = (4042.0 + 2263.0 * SQRT6) / 13824.0;
   static const double b13_7 = (-231278.0 + 40717.0 * SQRT6) / 69120.0;
   static const double b13_8 = (7947.0 - 2168.0 * SQRT6) / 11520.0;
   static const double b13_9 = (1048.0 - 542.0 * SQRT6) / 405.0;
   static const double b13_10 = (-1383.0 + 542.0 * SQRT6) / 720.0;
   static const double b13_11 = 2624.0 / 1053.0;
   static const double b13_12 = 3.0 / 1664.0;
   static const double b14_1 = -137.0 / 1296.0;
   static const double b14_6 = (5642.0 - 337.0 * SQRT6) / 864.0;
   static const double b14_7 = (5642.0 + 337.0 * SQRT6) / 864.0;
   static const double b14_8 = -299.0 / 48.0;
   static const double b14_9 = 184.0 / 81.0;
   static const double b14_10 = -44.0 / 9.0;
   static const double b14_11 = -5120.0 / 1053.0;
   static const double b14_12 = -11.0 / 468.0;
   static const double b14_13 = 16.0 / 9.0;
   static const double b15_1 = (33617.0 - 2168.0 * SQRT6) / 518400.0;
   static const double b15_6 = (-3846.0 + 31.0 * SQRT6) / 13824.0;
   static const double b15_7 = (155338.0 - 52807.0 * SQRT6) / 345600.0;
   static const double b15_8 = (-12537.0 + 2168.0 * SQRT6) / 57600.0;
   static const double b15_9 = (92.0 + 542.0 * SQRT6) / 2025.0;
   static const double b15_10 = (-1797.0 - 542.0 * SQRT6) / 3600.0;
   static const double b15_11 = 320.0 / 567.0;
   static const double b15_12 = -1.0 / 1920.0;
   static const double b15_13 = 4.0 / 105.0;
   static const double b16_1 = (-36487.0 - 30352.0 * SQRT6) / 279600.0;
   static const double b16_6 = (-29666.0 - 4499.0 * SQRT6) / 7456.0;
   static const double b16_7 = (2779182.0 - 615973.0 * SQRT6) / 186400.0;
   static const double b16_8 = (-94329.0 + 91056.0 * SQRT6) / 93200.0;
   static const double b16_9 = (-232192.0 + 121408.0 * SQRT6) / 17475.0;
   static const double b16_10 = (101226.0 - 22764.0 * SQRT6) / 5825.0;
   static const double b16_11 = - 169984.0 / 9087.0;
   static const double b16_12 = - 87.0 / 30290.0;
   static const double b16_13 =  492.0 / 1165.0;
   static const double b16_15 =  1260.0 / 233.0;

                           // e2 = 0, ..., e7 =0 //

   static const double e1 = -1911.0 / 109200.0;
   static const double e8 = 34398.0 / 109200.0;
   static const double e9 = -61152.0 / 109200.0;
   static const double e10 = 114660.0 / 109200.0;
   static const double e11 = -114688.0 / 109200.0;
   static const double e12 = -63.0 / 109200.0;
   static const double e13 = -13104.0 / 109200.0;
   static const double e14 = -3510.0 / 109200.0;
   static const double e15 = 39312.0 / 109200.0;
   static const double e16 = 6058.0 / 109200.0;

   double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16;
   double h12 = a2 * h;
   double h6 = a3 * h;

   k1 = (*f)(x0, *y);
   k2 = (*f)(x0+h12, *y + h12 * k1);
   k3 = (*f)(x0+a3*h, *y + h * ( b31*k1 + b32*k2) );
   k4 = (*f)(x0+a4*h, *y + h * ( b41*k1 + b43*k3) );
   k5 = (*f)(x0+a5*h, *y + h * ( b51*k1 + b53*k3 + b54*k4) );
   k6 = (*f)(x0+a6*h, *y + h * ( b61*k1 + b64*k4 + b65*k5) );
   k7 = (*f)(x0+a7*h, *y + h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6) );
   k8 = (*f)(x0+a8*h, *y + h * ( b81*k1 + b86*k6 + b87*k7) );
   k9 = (*f)(x0+a9*h, *y + h * ( b91*k1 + b96*k6 + b97*k7 + b98*k8) );
   k10 = (*f)(x0+a10*h, *y + h * ( b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8
                                                                + b10_9*k9 ) );
   k11 = (*f)(x0+a11*h, *y + h * ( b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9
                                                            + b11_10 * k10 ) );
   k12 = (*f)(x0+a12*h, *y + h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8
                                  + b12_9*k9 + b12_10 * k10 + b12_11 * k11 ) );
   k13 = (*f)(x0+a13*h, *y + h * ( b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8
                         + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12 ) );
   k14 = (*f)(x0+h, *y + h * ( b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8
            + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13 ) );
   k15 = (*f)(x0+a15*h, *y + h * ( b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8
            + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13 ) );
   k16 = (*f)(x0+h, *y + h * ( b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8
            + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13
                                                               + b16_15*k15) );
   *(y+1) = *y +  h * ( c1 * k1 + c8 * k8 + c9 * k9 + c10 * k10 + c11 * k11
                                         + c12 * k12 + c13 * k13 + c14 * k14 );
   return e1*k1 + e8*k8 + e9*k9 + e10*k10 + e11*k11 + e12*k12 + e13*k13
                                                 + e14*k14 + e15*k15 + e16*k16;
}

// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//   
//    //static const double err_exponent = 1.0 / 8.0;
// 
//    //double *scale;
//    //double *temp_y[2];
//    //double *err;
//    //double *yy;
//    //int *i;
//    //int *last_interval = 0;
//    
//    size_t nrows, ncols;
//    
//    mrows = mxGetM(phrs[0]);
//    ncols = mxGetN(phrs[0]);
//    
//    plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
//    
//    x = mxGetPr(phrs[0]);
//    y = mxGetPr(plhs[0]);
//    
//    
// }
