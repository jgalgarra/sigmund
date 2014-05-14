
function bino_new(dias)
    fn=mfilename();
    fnwp = mfilename('fullpath');
    fdat=strrep(fnwp, fn, '');
    Ma=dlmread(strcat(fdat,'..\input\1dat_new_a.txt'));
    Mb=dlmread(strcat(fdat,'..\input\1dat_new_b.txt'));

    r1=Ma(4);  % nacimientos año por bicho
    r1d=Ma(5);
    period=365
    N1_0=Ma(2);
    N1log=N1_0;
    K1=Ma(3); % capcidad limite
    %beta12=Ma(1)*K1/r1
    beta12=Ma(1)
    
    %beta12period=beta12/period;

    r2=Mb(4);  % nacimientos año por bicho
    r2d=Mb(5);
    N2_0=Mb(2);
    N2log=N2_0;
    K2=Mb(3); % capcidad limite
    %beta21=Mb(1)*K2/r2
    
    beta21=Mb(1)


    for k=1:dias-1
        
        if k==366
            z=N1logAnt
            w=N2logAnt
        end
        
        N1logAnt = N1log;
        N2logAnt = N2log;
        rmut1 = beta12*N2logAnt;
        r1b = r1+rmut1;
        if (r1b==0)
            p1period = 0;
            incN1log = 0;
            t1b=0;
        else
            t1b = -1*abs((r1d-r1b))*(N1logAnt/K1)+r1b;
            p1period = rp(t1b,period);
            incN1log=binornd(N1logAnt,abs(p1period));
        end
        if (r1d==0)
            decN1log = 0;
        else
            p1dperiod = rp(r1d,period);
            decN1log = binornd(N1logAnt,abs(p1dperiod));
        end
        N1log=N1logAnt+sign(p1period)*incN1log-decN1log;
        
        rmut2= beta21*N1logAnt;
        r2b = r2+rmut2;
        if (r2b==0)
            p2period = 0;
            incN2log = 0;
            t2b=0;
        else
            t2b = -1*abs(r2d-r2b)*(N2logAnt/K2)+r2b;
            p2period = rp(t2b,period);
            incN2log=binornd(N2logAnt,abs(p2period));
        end
        if (r2d==0)
            decN2log = 0;
        else
            p2dperiod = rp(r2d,period);
            decN2log = binornd(N2logAnt,abs(p2dperiod));
        end
        N2log=N2logAnt+sign(p2period)*incN2log-decN2log;
        
        N1(k) = N1log;
        N2(k) = N2log;
        r1eq(k)=t1b-r1d;
        r2eq(k)=t2b-r2d;

        
    end

    hold on;
        subplot(2, 1, 1)
        plot(N1,'-r','LineWidth',2);
        hold on;
        plot(N2,'-k','LineWidth',2);
        title(sprintf('N1inic=%d K1=%d r1b=%0.2f r1d=%0.2f beta12=%0.6f (red); N2inic=%d K2=%d r2b=%0.2f r2d=%0.2f beta21=%0.8f (black)',N1_0,K1,r1,r1d,beta12,N2_0,K2,r2,r2d,beta21));
        hold on;
        subplot(2, 1, 2)
        %plot(k,t1b-r1d,'-r');
        %hold on;
        plot(r1eq,'--r','LineWidth',1);
        hold on;
        plot(r2eq,'--k','LineWidth',1);
        title(sprintf('r1eq (red); r2eq (black)'));
        
