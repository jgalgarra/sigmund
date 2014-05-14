%% HOW DO I DO THAT IN MATLAB SERIES?
% In this series, I am answering questions that students have asked
% me about MATLAB.  Most of the questions relate to a mathematical 
% procedure.

%% TOPIC
% How do I solve a nonlinear equation?

%% SUMMARY

% Language : Matlab 2008a; 
% Authors : Autar Kaw; 
% Mfile available at
% http://numericalmethods.eng.usf.edu/blog/integration.m; 
% Last Revised : March 28, 2009; 
% Abstract: This program shows you how to solve a nonlinear equation.
clc
clear all

%% INTRODUCTION

disp('ABSTRACT')
disp('   This program shows you how to solve')
disp('   a nonlinear equation')
disp(' ')
disp('AUTHOR')
disp('   Autar K Kaw of http://autarkaw.wordpress.com')
disp(' ')
disp('MFILE SOURCE')
disp('   http://numericalmethods.eng.usf.edu/blog/nonlinearequation.m')
disp(' ')
disp('LAST REVISED')
disp('   April 11, 2009')
disp(' ')

%% INPUTS
% Solve the nonlinear equation x^3-15*x^2+47*x-33=0
% Define x as a symbol
syms x
% Assigning the left hand side of the equation f(x)=0
f=x^3-15*x^2+47*x-33;
%% DISPLAYING INPUTS

disp('INPUTS')
func=['  The equation to be solved is ' char(f), '=0'];
disp(func)
disp('  ')

%% THE CODE

% Finding the solution of the nonlinear equation
soln=solve(f,x);
solnvalue=double(soln);

%% DISPLAYING OUTPUTS

disp('OUTPUTS')
for i=1:1:length(solnvalue)
fprintf('\nThe solution# %g is %g',i,solnvalue(i))
end
disp('  ')