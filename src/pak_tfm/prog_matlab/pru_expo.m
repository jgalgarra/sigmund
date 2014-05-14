
function pru_expo(dias,N0,ranual)
    
    rdiaria = rp(ranual,365)
    N1log=N0;
    N1bin=N0;
    for k=1:dias
        
        N1logAnt = N1log;
        N1binAnt = N1bin;
        incN1bin=binornd(N1binAnt,abs(rdiaria));
        N1bin = N1binAnt + sign(rdiaria)*incN1bin;
        incN1log=N1logAnt*exp(rdiaria);
        N1log = incN1log;
        N1(k)=N1log;
        N1b(k)=N1bin;
    end
    
    N1log
    N1bin
    hold on;
    plot(N1,'-r','LineWidth',2);
    hold on;
    plot(N1b,'-b','LineWidth',2);
     
        
