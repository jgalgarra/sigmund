
function bino_May(dias)
    fn=mfilename();
    fnwp = mfilename('fullpath');
    fdat=strrep(fnwp, fn, '');
    Ma=dlmread(strcat(fdat,'..\input\1dat_a.txt'));
    Mb=dlmread(strcat(fdat,'..\input\1dat_b.txt'));

    r1=Ma(4);  % nacimientos año por bicho
    period=365
    r1period=r1/period;
    N1_0=Ma(2);
    N1log=N1_0;
    K1=Ma(3); % capcidad limite
    %beta12=Ma(1)*K1/r1
    beta12=Ma(1)
    
    %beta12period=beta12/period;

    r2=Mb(4);  % nacimientos año por bicho
    r2period=r2/period;
    N2_0=Mb(2);
    N2log=N2_0;
    K2=Mb(3); % capcidad limite
    %beta21=Mb(1)*K2/r2
    
    beta21=Mb(1)
    beta21period=beta21/period;

    for k=1:dias
        %incN1logexp=binornd(N1log,1-exp(-1*r1period));
        %incN1log=binornd(round(N1log^2/K1),1-exp(-1*r1period));
        %incN1log_esp2=binornd(round(beta12*N1log*N2log/K1),1-exp(-1*r1period));
        %N1log = N1log + incN1logexp - incN1log + incN1log_esp2;
        
        %term1 = N1log - N1log^2/K1 + beta12*N1log*N2log/K1;
        %incN1log=binornd(round(abs(term1)),1-exp(-1*r1period));
        %N1log=N1log+sign(incN1log)*incN1log;

        term1 = 1 - N1log/K1 + beta12*N2log/K1;
        p1period = rp(abs(r1*term1),period);
        incN1log=binornd(round(N1log),1-exp(-1*p1period));
        N1log=N1log+sign(incN1log)*incN1log;

        
        %incN2logexp=binornd(N2log,1-exp(-1*r2period));
        %incN2log=binornd(round(N2log^2/K2),1-exp(-1*r2period));
        %incN2log_esp2=binornd(round(beta21*N2log*N1log/K2),1-exp(-1*r2period));
        %N2log = N2log + incN2logexp - incN2log + incN2log_esp2;
        
        %term2 = N2log - N2log^2/K2 + beta21*N2log*N1log/K2;
        %incN2log=binornd(round(abs(term2)),1-exp(-1*r2period));
        %N2log=N2log+sign(incN2log)*incN2log;
        
        term2 = 1 - N2log/K2 + beta21*N1log/K2;
        p2period = rp(abs(r2*term2),period);
        incN2log=binornd(round(N2log),1-exp(-1*p2period));
        N2log=N2log+sign(incN2log)*incN2log;

        plot(k,N1log,'-.g');
        hold on;
        plot(k,N2log,'-.k');
    end

    title(sprintf('N1inic=%d K1=%d r1=%0.2f beta12=%0.2f (red); N2inic=%d K2=%d r2=%0.2f beta21=%0.2f (black)',N1_0,K1,r1,beta12,N2_0,K2,r2,beta21));

