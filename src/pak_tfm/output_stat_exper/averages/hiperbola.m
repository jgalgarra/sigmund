x = linspace(250,4500);
y=  (x.^(-1.482))
c = exp(16.58)
y = c*y
plot(x,y,'.k');
hold on;

%xlimit = ([0 5000]);
%ylimit = ([0 5000]);