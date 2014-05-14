Min1=dlmread('probslice_nexper_1x1_Binary_Logistic_abs_xcentral655_ycentral1030_numpoints40_nexper_300_years_100.txt');
Min1(:,1)=int32(Min1(:,1));
Min1(:,2)=int32(Min1(:,2));
ycoord = Min1(1,2,1);
xsup = max(Min1(:,1));
xinf = min(Min1(:,1));
xmax=800;
xmin =400;
for j=xsup+1:xmax
    fila = [j,ycoord,1.0];
    Min1 = vertcat(Min1,fila);
end

for j=xmin:xinf-1
    fila = [j,ycoord,0.0];
    Min1 = vertcat(fila,Min1);
end
figure
plot(Min1(:,1),Min1(:,3),'-r');
ylim = [0, 2];
xlim = [500, 1100];

hold on;

