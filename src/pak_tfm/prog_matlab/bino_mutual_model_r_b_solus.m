%  Mutualist model modified with r-b verhulst formulation

clear;
r1=-0.005;  % nacimientos año por bicho
r1month=r1/12;
N1_0=1100;
N1log_ant=N1_0;
N1log=N1_0;
K1=10000; % capcidad limite en vacío
coef_b1_month = r1/(12*K1);
b12=0.00001850;
b12month=b12/12;
coef_roz1 = 0.0000055/12;
incroz12 = 0.0000000010/12;

r2=-0.005;  % nacimientos año por bicho
r2month=r2/12;
N2_0=800;
N2log_ant=N2_0;
N2log = N2_0;
K2=8000; % capcidad limite en vacío
b21=0.00001750;
b21month=b21/12;
coef_b2_month = r2/(12*K2);
coef_roz2 = 0.00001/12;
incroz21 = 0.0000000010/12;

constroz1 = 0.00004;
constroz2 = 0.00018;


for x=0:5000
    y = (r1month + b12month * x)/(coef_roz1 + incroz12 * b12month * x);
    plot(x,y)
    hold on
end
for y=0:5000
    x = (r2month + b21month * y)/(coef_roz2 + incroz21 * b21month * y);
    hold on
end