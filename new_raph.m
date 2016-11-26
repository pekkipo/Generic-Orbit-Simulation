function []= new_raph(xold,n)
syms x
fx1=f(xold);
f_x=diff(f(x));
fdd=subs(f_x,xold);
for i=1:n
    xnew = xold - (fx1/fdd);
    fx1=f(xnew);
    fdd=subs(f_x,xnew);
    err=abs((xnew-xold)/xnew)*100;
    xold=xnew;
    disp('the new value of x is');disp(xnew)
    disp('the relative error in percentage is');disp(err)
end
end
    
