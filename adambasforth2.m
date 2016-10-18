function [t, y] = adambasforth2( fmodel, range, initial_y, numSteps )
% [ x, y ] = ab2 ( f_ode, xRange, yInitial, numSteps ) uses
% Adams-Bashforth second-order method to solve a system 
% of first-order ODEs dy/dx=f_ode(x,y).
% f = name of an m-file with signature
%    fValue = f_ode(x,y)
% to compute the right side of the ODE as a column vector
% force_model m-file in my case
%
% range = [x1,x2] where the solution is sought on x1<=x<=x2
% yInitial = column vector of initial values for y at x1
% numSteps = number of equally-sized steps to take from x1 to x2
% x = row vector of values of x
% y = matrix whose k-th row is the approximate solution at x(k).

x(1) = xRange(1);
h = ( xRange(2) - xRange(1) ) / numSteps;
y(:,1) = yInitial;

k = 1;
  fValue =  f_ode( x(k), y(:,k) );
  xhalf = x(k) + 0.5 * h;
  yhalf = y(:,k) + 0.5 * h * fValue;
  fValuehalf = f_ode( xhalf, yhalf );

  x(1,k+1) = x(1,k) + h;
  y(:,k+1) = y(:,k) + h * fValuehalf;

for k = 2 : numSteps
  fValueold=fValue;
  fValue = f_ode( x(k), y(:,k) );
  x(1,k+1) = x(1,k) + h;
  y(:,k+1) = y(:,k) + h * ( 3 * fValue - fValueold ) / 2;
end
end

