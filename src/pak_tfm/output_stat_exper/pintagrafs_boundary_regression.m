Ma= dlmread('bd_nexper_1x1_00005_Logistic_abs_Pl_r_b_-0.050000_0.000050_Pol_r_b_-0.090000_0.000110.txt');
hold on;
lonserie = size(Ma);
limit = 1+(lonserie(1)-1)/2;
h = plot(Ma(1:limit,1),Ma(1:limit,2),'.r');
xlim([0 max(Ma(:,1))])
ylim([0 max(Ma(:,2))])
set(h, 'Markersize',12);

h = plot(Ma(limit:lonserie(1),1),Ma(limit:lonserie(1),2),'.b');
xlim([0 max(Ma(:,1))])
ylim([0 max(Ma(:,2))])
set(h, 'Markersize',12);

%{
Ma=dlmread('bd_nexper_1x1_Logistic_u_Pl_r_b_-0.050000_0.000080_Pol_r_b_-0.090000_0.000110.txt');
hold on;
h = plot(Ma(:,1),Ma(:,2),'.b');
xlim([0 max(Ma(:,1))])
ylim([0 max(Ma(:,2))])
set(h, 'Markersize',12);
%}

%{
Mbichos=dlmread('popevo_nexper_1x1_normal.txt');
hold on;
h = plot(Mbichos(:,1),Mbichos(:,2),'.r');
set(h, 'Markersize',1);


Mbichos=dlmread('popevo_chunga1x1_survival.txt');
hold on;
h = plot(Mbichos(:,1),Mbichos(:,2),'.g');
set(h, 'Markersize',1);

Mbichos=dlmread('popevo_chunga1x1_death.txt');
hold on;
h = plot(Mbichos(:,1),Mbichos(:,2),'.b');
set(h, 'Markersize',1);
%}
%{
hold on;
Mb=dlmread('bd_nexper_1x1_00010_Logistic_abs_Pl_r_b_-0.050000_0.000100_Pol_r_b_-0.090000_0.000110.txt');
plot(Mb(:,1),Mb(:,2),'.r');



hold on;
Mc=dlmread('bd_nexper_1x1_00005_Logistic_abs_Pl_r_b_-0.050000_0.000050_Pol_r_b_-0.090000_0.000110.txt');%
plot(Mc(:,1),Mc(:,2),'.g');


hold on;
Md=dlmread('bd_nexper_1x1_b00009_Logistic_abs_Pl_r_b_-0.050000_0.000080_Pol_r_b_-0.090000_0.000090.txt');%
plot(Md(:,1),Md(:,2),'.g');
%}