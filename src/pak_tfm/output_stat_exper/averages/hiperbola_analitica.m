% a =0.02

x = linspace(200,4500,400);
exponente=-1.2939;
y=  (x.^(exponente));
c =  989.0388/692.6974^exponente;
y = c*y;
plot(x,y,'-b');
hold on;
return
%xlimit = ([0 5000]);
%ylimit = ([0 5000]);

% a = 0.015
exponente=-1.21;
y=  (x.^(exponente));
c = 1381.7/906.5^exponente;
y = c*y;
plot(x,y,'.r');
hold on;
