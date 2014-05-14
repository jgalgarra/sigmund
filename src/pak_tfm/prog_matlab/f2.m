function [f] = f2(x)
% f2(x) = 0, with x = [x(1);x(2)]
% represents a system of 2 non-linear equations
f1 = r(1)+b(1)*x(2)-alpha(1)*x(1)-c(1)*b(1)*x(1)*x(2);
f2 = r(2)+b(2)*x(1)-alpha(2)*x(2)-c(2)*b(2)*x(2)*x(1);
f = [f1;f2];
% end function

