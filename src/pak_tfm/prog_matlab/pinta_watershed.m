% Genera y pinta una watershed para dos especies 
%
% N2ini = (N2crit * N1ini^-B/A) / N1crit^-B/A
%
% Guarda los valores de N1 y N2 en el fichero watershed.txt


% Los valores iniciales han sido calculados previamente con el
% script calK1K2.m

N1crit = 692;
N2crit = 989;
limsuppop = 5000;
liminfpop = 100;
exponente = -1.2944 

N1 = liminfpop:10:limsuppop;
exponente = -B/A;
N2 = (N2crit/N1crit^exponente)*N1.^exponente;
ndatos = length(N1);
for i = 1:ndatos
    results(i,:) = [N1(i),N2(i)];
end 
plot(N1,N2);
axis([0 5000 0 5000]);
dlmwrite('watershed.txt',results,'newline','pc');