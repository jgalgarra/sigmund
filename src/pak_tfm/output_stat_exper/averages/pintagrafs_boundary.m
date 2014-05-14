
M1=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM1.txt');
M2=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM2.txt');
M3=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM3.txt');
M4=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM4.txt');
M5=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM5.txt');
M6=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM6.txt');
M7=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM7.txt');
M8=dlmread('bd_vexper_1x1_Verhulst_Binary_1.0_0.1_NUM8.txt');

%{
mediax = round(1/8* ((M1(:,1)+M2(:,1)+M3(:,1)+M4(:,1)+M5(:,1)+M6(:,1)+M7(:,1)+M8(:,1)+M9(:,1)+M10(:,1))));
mediay = round(1/8* ((M1(:,2)+M2(:,2)+M3(:,2)+M4(:,2)+M5(:,2)+M6(:,2)+M7(:,2)+M8(:,2)+M9(:,2)+M10(:,2))));
%}
mediax = round(1/8* ((M1(:,1)+M2(:,1)+M3(:,1)+M4(:,1)+M5(:,1)+M6(:,1)+M7(:,1)+M8(:,1))));
mediay = round(1/8* ((M1(:,2)+M2(:,2)+M3(:,2)+M4(:,2)+M5(:,2)+M6(:,2)+M7(:,2)+M8(:,2))));



hold on;
h = plot(mediax,mediay,'.r');
%%set(h,'LineWidth',2);

mediax
mediay

return

Mb1=dlmread('bd_vexper_1x1_r5_Verhulst_Binary_1.0_0.1_NUM1.txt');
Mb2=dlmread('bd_vexper_1x1_r5_Verhulst_Binary_1.0_0.1_NUM0.txt');
%{
Mb3=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM3.txt');
Mb4=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM4.txt');
Mb5=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM5.txt');
Mb6=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM6.txt');
Mb7=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM7.txt');
Mb8=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM8.txt');
Mb9=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM9.txt');
Mb10=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1_NUM9.txt');
%}

mediax = round(1/2*((Mb1(:,1)+Mb2(:,1))) );
mediay = round(1/2*((Mb1(:,2)+Mb2(:,2))) );

hold on;
t = plot(mediax,mediay,'-b');
return

Mc1=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM1.txt');
Mc2=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM2.txt');
Mc3=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM3.txt');
Mc4=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM4.txt');
Mc5=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM5.txt');

mediax = round(1/5*(Mc1(:,1)+ Mc2(:,1)+ Mc3(:,1)+ Mc4(:,1)+ Mc5(:,1)));
mediay = round(1/5*(Mc1(:,2)+ Mc2(:,2)+ Mc3(:,2)+ Mc4(:,2)+ Mc5(:,2)));

hold on;
q = plot(mediax,mediay,'-k');


Md1=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.96_0.1_NUM2.txt');
Md2=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM2.txt');
Md3=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.98_0.1_NUM3.txt');

mediax = round(1/3*(Md1(:,1)+ Md2(:,1)+ Md3(:,1)));
mediay = round(1/3*(Md1(:,2)+ Md2(:,2)+ Md3(:,2)));

hold on;
q = plot(mediax,mediay,'-g');

Me1=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.97_0.1_NUM1.txt');
Me2=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.97_0.1_NUM2.txt');
Me3=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.97_0.1_NUM3.txt');

mediax = round(1/3*(Me1(:,1)+ Me2(:,1)+ Me3(:,1)));
mediay = round(1/3*(Me1(:,2)+ Me2(:,2)+ Me3(:,2)));


hold on;
r = plot(mediax,mediay,'--m');
%%xlim([0 max(Ma(:,1))])
%%ylim([0 max(Ma(:,2))])
%%set(h, 'Markersize',12);

xlim([0 1400]);
ylim([0 1300]);

%{
Mb=dlmread('bd_nexper_1x1_Logistic_abs_Pl_r_b_-0.050000_0.000080_Pol_r_b_-0.090000_0.000110_NUM1.txt');
hold on;
h = plot(Mb(:,1),Mb(:,2),'-r');

Mc=dlmread('bd_nexper_1x1_Logistic_abs_Binary_0.99_0.1.txt');
hold on;
h = plot(Mc(:,1),Mc(:,2),'-b');
%}

%{
M10=dlmread('bd_nexper_1x1_Logistic_abs_Pl_r_b_-0.050000_0.000080_Pol_r_b_-0.090000_0.000110_NUM10.txt');
hold on;
h = plot(M10(:,1),M10(:,2),'-b');
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