Ma=dlmread('z_nexper_1x1_GaussianLogistic_abs_sd0.05_ini80.0_fin100.0_numsims200_pasos_20.txt');
plot(Ma(:,1),Ma(:,2),'--r');
hold on;
Mb=dlmread('z_nexper_1x1_GaussianLogistic_abs_sd0.1_ini80.0_fin100.0_numsims200_pasos_20.txt');
plot(Ma(:,1),Mb(:,2),'--b');
hold on;
Mc=dlmread('z_nexper_1x1_GaussianLogistic_abs_sd0.2_ini80.0_fin100.0_numsims200_pasos_20.txt');%
plot(Ma(:,1),Mc(:,2),'--g');