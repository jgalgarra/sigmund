function ydot=diff_May(t,y,r1,K1,b12,r2,K2,b21)
       %     right hand side of logistic equation for a matlab numerical
       %     solution.

       
       ydot=[r1*y(1)*(1-(y(1)/K1)+b12*y(2)/K1);r2*y(2)*(1-(y(2)/K2)+b21*y(1)/K2)];