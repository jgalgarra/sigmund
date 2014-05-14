function [x,iter]=newton(x0,f,fp)
% newton-raphson algorithm
N = 100; eps = 1.e-5; % define max. no. iterations and error
maxval = 10000.0; % define value for divergence
xx = x0;
while (N>0)
    xn = xx-f(xx)/fp(xx);
    if abs(f(xn))<eps
        x=xn;iter=100-N;
        return;
    end;
    if abs(f(xx))>maxval
        disp(['iterations = ',num2str(iter)]);
        error('Solution diverges');
        break;
    end;
    N = N - 1;
    xx = xn;
end;
error('No convergence');
return;
% end function