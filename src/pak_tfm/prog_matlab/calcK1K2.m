% Este script genera una gráfica con las soluciones del modelo para dos
% especies

b12 = 0.000041850;

c1 = 0.00004;
alpha1 = 0.000035;
%r1 =0.02-0.036;
r1=0.02-0.036;

b21=0.00008750; % dos raices
%b21= 0.000063185; % una raiz
%b21 = 0.00004; % raices complejas
c2 = 0.0001;
alpha2 = 0.000035;
r2 =0.05-0.07;

K1_A = (c2*b21*alpha1+c1*b12*b21);
K1_B = (alpha1*alpha2 + c1*b12*r2 - c2*b21*r1 -b12*b21);
K1_C = -1*(r1*alpha2 + b12* r2);

raiz_K1_1 = (-1*K1_B+sqrt(K1_B^2-4*K1_A*K1_C))/(2*K1_A)
raiz_K1_2 = (-1*K1_B-sqrt(K1_B^2-4*K1_A*K1_C))/(2*K1_A)


K2_A = (c1*b12*alpha2+c2*b21*b12);
K2_B = (alpha2*alpha1 + c2*b21*r1 - c1*b12*r2 -b21*b12);
K2_C = -1*(r2*alpha1 + b21* r1);


raiz_K2_1 = (-1*K2_B+sqrt(K2_B^2-4*K2_A*K2_C))/(2*K2_A)
raiz_K2_2 = (-1*K2_B-sqrt(K2_B^2-4*K2_A*K2_C))/(2*K2_A)
 
N1 = linspace(0,5000);
N2=  (alpha1*N1-r1)./(b12 - c1*b12*N1);

plot(N1,N2,'-r');
hold on;
xlimit = ([0 5000]);
ylimit = ([0 5000]);

N2 = linspace(0,5500);
N1=  (alpha2*N2-r2)./((b21 - c2*b21*N2));

plot(N1,N2,'-b');
hold on;

% Calculo de b21 critico

m = c2*alpha1+c1*b12;
n = c2*r1+b12;
p = alpha1*alpha2 + c1*b12*r2;
c = K1_C;

BA = n^2;
BB = -4*m*c-2*n*p;
BC = p^2;


raiz_b21_1 = (-1*BB+sqrt(BB^2-4*BA*BC))/(2*BA)
raiz_b21_2 = (-1*BB-sqrt(BB^2-4*BA*BC))/(2*BA)

% Calculo de condicion de equilibrio
N1crit=raiz_K1_2
N2crit=raiz_K2_2
%N1ini = 553 %548
%N2ini = 1323%1385
N1ini = 553
N2ini = 1323
A=((r1+b12*N2ini) - N1ini*(alpha1+c1*b12*N2ini))
B=-1*((r2+b21*N1ini) - N2ini*(alpha2+c2*b21*N1ini))

%N1e=N1ini+((N1crit-N1ini)*(1-1/exp(1)))
%kanual=A*(1-exp(-1)) / log(N1e/N1ini)
%kdia =  exp(kanual/365)-1;


%syms x
% Assigning the left hand side of the equation f(x)=0
G=1000;
f=@(x) G*(log(N1crit/N1ini)* ((r2+b21*N1ini) - x *(alpha2+c2*b21*N1ini)) - log(N2crit/x) * ((r1+b12*x) - N1ini*(alpha1+c1*b12*x)));
% Numerical solution
h = fsolve(f,1500)


