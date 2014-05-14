bp(1,1) = 0.000021850;
bp(1,2) = 0.000051;
bp(2,1) = 0.000001;
bp(2,2) = 0;

cp(1) = 0.00004;
cp(2) = 0.00003;
alphap(1) = 0.000035;
alphap(2) = 0.000024;
rp(1)=0.02-0.036;
%rp(2)=0.02-0.021;
rp(2)=0;

ba(1,1)= 0.0000265;
ba(1,2)= 0.00004;
ba(2,1)= 0.00002950;
ba(2,2)= 0;
ca(1) = 0.0001;
ca(2) = 1000;
alphaa(1)=0.000035;
alphaa(2)=1000;
ra(1)=0.05-0.07;
ra(2)=0;

B = [-alphap(1) 0 bp(1,1) bp(2,1);
     0 -alphap(2) bp(1,2) bp(2,2);
     ba(1,1) ba(2,1) -alphaa(1) 0;
     ba(1,2) ba(2,2) 0 -alphaa(2)];
m = [-rp(1);-rp(2);-ra(1);-ra(2)]; 
gxmin = linsolve(B,m)

G=10;

x0=gxmin;
inda = 2;
f2 =@(x) [rp(1)+bp(1,1)*x(1+inda)+bp(2,1)*x(2+inda)-x(1)*(alphap(1)+cp(1)*(bp(1,1)*x(1+inda)+bp(2,1)*x(2+inda)));
          rp(2)+bp(1,2)*x(1+inda)+bp(2,2)*x(2+inda)-x(2)*(alphap(2)+cp(2)*(bp(1,2)*x(1+inda)+bp(2,2)*x(2+inda)));
          ra(1)+ba(1,1)*x(1)+ba(2,1)*x(2)-x(1+inda)*(alphaa(1)+ca(1)*(ba(1,1)*x(1)+ba(2,1)*x(2)));
          ra(2)+ba(1,2)*x(1)+ba(2,2)*x(2)-x(2+inda)*(alphaa(2)+ca(2)*(ba(1,2)*x(1)+ba(2,2)*x(2)))];
f = @(x)[G*f2(x)];
xmin = fsolve(f,x0)

x0=[1/cp(1);1/cp(2);1/ca(1);1/ca(2)];
xmax = fsolve(f,x0)