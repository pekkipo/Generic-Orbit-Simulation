
function [y, t] = adambashforth8(f,vrange,y0,n, step)

    %h = length(vrange)  / n;
    h = step;
    h120960 = h /120960;

    y(:,1) = y0;
    t(1) = vrange(1);

    m = min(7,n);

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
        %y(:,i+1) = y(:,i) + (55 * s(:,i) - 59 * s(:,i-1) + 37 * s(:,i-2) - 9 * s(:,i-3)) * h24; % predictor
        y(:,i+1) = y(:,i) + (434241 * s(:,i) - 1152169 * s(:,i-1) + 2183877 * s(:,i-2) - 2664477 * s(:,i-3) + 2102243 * s(:,i-4) - 1041723 * s(:,i-5) + 295767 * s(:,i-6) - 36799 * s(:,i-7)) * h120960; % predictor
        t(i+1) = t(i) + h;
       % y(:,i+1) = y(:,i) + (9 * f(t(i+1), y(:,i+1)) + 19 * s(:,i) - 5 * s(:,i-1) + s(:,i-2)) * h24; % corrector
        y(:,i+1) = y(:,i) + (36799 * f(t(i+1), y(:,i+1)) + 139849 * s(:,i) - 121797 * s(:,i-1) + 123133*s(:,i-2) - 88547*s(:,i-3) + 41499*s(:,i-4) - 11351*s(:,i-5) + 1375*s(:,i-6)) * h120960; % corrector
    end;
end
