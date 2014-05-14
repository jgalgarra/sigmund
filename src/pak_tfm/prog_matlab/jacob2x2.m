function [J] = jacob2x2(x)
% Evaluates the Jacobian of a 2x2
% system of non-linear equations
J(1,1) = -alpha(1)-c(1)*b(1)*x(2); J(1,2) = b(1)-c(1)*b(1)*x(1);
J(2,1) = b(2)-c(2)*b(2)*x(2); J(2,2) = -alpha(2)-c(2)*b(2)*x(1);
% end function

