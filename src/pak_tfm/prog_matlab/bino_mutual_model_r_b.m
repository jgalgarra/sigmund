%  Mutualist model modified with r-b verhulst formulation

clear;
r1=-0.005;  % nacimientos a�o por bicho
r1month=r1/12;
N1_0=1100;
N1log_ant=N1_0;
N1log=N1_0;
K1=10000; % capcidad limite en vac�o
coef_b1_month = r1/(12*K1);
b12=0.00001850;
b12month=b12/12;
coef_roz1 = 0.0000055/12;
incroz12 = 0.0000000010/12;

r2=-0.005;  % nacimientos a�o por bicho
r2month=r2/12;
N2_0=800;
N2log_ant=N2_0;
N2log = N2_0;
K2=8000; % capcidad limite en vac�o
b21=0.00001750;
b21month=b21/12;
coef_b2_month = r2/(12*K2);
coef_roz2 = 0.00001/12;
incroz21 = 0.0000000010/12;

constroz1 = 0.00004;
constroz2 = 0.00018;

meses=4000;
for k=1:meses
    %r1eq = r1month+r1month*b12month*N2log/K1-r1month*N1log/K1;
    termeq1 = b12month*N2log_ant;
    r1eq = r1month+termeq1;
    
    %r1eq =r1month;
    
    %roz1 = (coef_roz1+incroz12*N2log)*N1log_ant;
    roz1 = (coef_roz1+termeq1*constroz1)*N1log_ant;
    rspneq1 = r1eq - roz1;
    incN1=binornd(N1log_ant,1-exp(-1*abs(rspneq1)));
    %decN1 = binornd(N1log,1-exp(-1*abs(roz1)));
    N1log = N1log_ant +sign(rspneq1)*incN1;
    K1cal = r1eq/(coef_roz1+termeq1*constroz1);
   
    %r2eq = r2month+r2month*b21month*N1log/K2-r2month*N2log/K2;
    %r2eq = r2month+b21month*N1log_ant;
    termeq2 = b21month*N1log_ant;
    r2eq = r2month+termeq2;
    
    %r2eq = r2month;
    
    %roz2 = (coef_roz2+incroz21*N1log_ant)*N2log;
    roz2 = (coef_roz2+termeq2*constroz2)*N2log_ant;
    rspneq2 = r2eq - roz2;
    incN2=binornd(N2log_ant,1-exp(-1*abs(rspneq2)));
    %decN2=binornd(N2log,1-exp(-1*abs(roz2)));
    N2log = N2log_ant +sign(rspneq2)*incN2;
    K2cal = r2eq/(coef_roz2+termeq2*constroz2);
    
    N1log_ant = N1log;
    N2log_ant =N2log;
    hold on;
    plot(k,N1log,'--b','LineWidth',1);
    %plot(k,rspneq1,'-.r','LineWidth',4);
    %plot(k,roz1,'-.k','LineWidth',4);
    plot(k,K1cal,'-r','LineWidth',1);
    %scatter(N1log,N2log,'.');
    
    
    hold on;
    %plot(k,roz2,'-b','LineWidth',4);
    %plot(k,rspneq2,'-b','LineWidth',4);
    %plot(k,N2log,'--b','LineWidth',4);
    %plot(k,K2cal,'-b','LineWidth',1);
   % plot(k,incN2,'-g','LineWidth',4);


end
    
%title(sprintf('N1inic=%d K1=%d r1=%0.2f b12=%0.2f (red); N2inic=%d K2=%d r2=%0.2f b21=%0.2f (blue)',N1_0,K1,r1,b12,N2_0,K2,r2,b21));

