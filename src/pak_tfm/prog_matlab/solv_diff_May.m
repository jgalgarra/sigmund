tspan=[0 1500];
       y0=[1000 2000];
       r1=0.1/12;
       K1=10000;
       b12=0.1/12;
       
       r2=0.1/12;
       K2=8000;
       b21=0.1/12;
       b23=0.00;
       K3=100;
       b31=0.0;
       b32=0;
       r3=0.0;
       options = odeset('reltol',1e-6);
       tic;
       [t,y] = ode23(@mutual, tspan,y0,options,r1,K1,b12,r2,K2,b21);
       toc;
       plot(t,y,'-.k','LineWidth',2), xlabel('Time'), ylabel('Population')