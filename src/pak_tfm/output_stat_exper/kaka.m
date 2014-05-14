[X,Y] = meshgrid(-8:.5:8);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;

% surface in 3D
figure;
surf(Z,'EdgeColor','None');