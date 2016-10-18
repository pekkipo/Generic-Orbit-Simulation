% Usage: [y t] = abm4(f,a,b,ya,n) or y = abm4(f,a,b,ya,n)
% Adams-Bashforth-Moulton 4-th order predictor-corrector method for initial value problems
% It uses 
% Adams-Bashforth 4-step method as a precdictor,
% Adams-Moulton 3-step method as a corrector, and
% Runge-Kutta method of order 4 as a starter
%
% Input:
% f - Matlab inline function f(t,y)
% a,b - interval
% y0 - initial condition
% n - number of subintervals (panels)
%
% Output:
% y - computed solution
% t - time steps
%
% Examples:
% [y t]=abm4(@myfunc,0,1,1,10);          here 'myfunc' is a user-defined function in M-file
% y=abm4(inline('sin(y*t)','t','y'),0,1,1,10);
% f=inline('sin(y(1))-cos(y(2))','t','y');
% y=abm4(f,0,1,1,10);

function [y, t] = adambashforth4(f,vrange,y0,n, step)

    %h = length(vrange)  / n;
    h = step;
    h24 = h / 24;

    y(:,1) = y0;
    t(1) = vrange(1);

    m = min(3,n);

    for i = 1 : m % start-up phase, using Runge-Kutta of order 4
        t(i+1) = t(i) + h;
        s(:,i) = f(t(i), y(:,i));
        s2 = f(t(i) + h / 2, y(:,i) + s(:,i) * h /2);
        s3 = f(t(i) + h / 2, y(:,i) + s2 * h /2);
        s4 = f(t(i+1), y(:,i) + s3 * h);
        y(:,i+1) = y(:,i) + (s(:,i) + s2+s2 + s3+s3 + s4) * h / 6;
    end;

    for i = m + 1 : n % main phase
        s(:,i) = f(t(i), y(:,i));
        y(:,i+1) = y(:,i) + (55 * s(:,i) - 59 * s(:,i-1) + 37 * s(:,i-2) - 9 * s(:,i-3)) * h24; % predictor
        t(i+1) = t(i) + h;
        y(:,i+1) = y(:,i) + (9 * f(t(i+1), y(:,i+1)) + 19 * s(:,i) - 5 * s(:,i-1) + s(:,i-2)) * h24; % corrector
    end;
end