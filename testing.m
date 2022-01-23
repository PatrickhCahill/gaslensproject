%% BOLTZMANN LENS
%Sets up symbolic function environment
syms a b

%Constant and other parameters %% Note that initial ray conditions are in
%trace() function
alpha = 1.642e-30; %polarizability
eps0 = 8.854187817e-12; 
k = 1.38064852e-23; %boltzmann
avagadro =6.02214086e23 ;

t = 0.1;
gasmass = 6000;
temp0 = 90;
molarmass = 39.948;
N = (gasmass/molarmass)*avagadro;
particlemass = molarmass/avagadro/1000;

beta0 = particlemass/(2*k*temp0);

%Input pressure function and symbolic computes refractive index and other
%functions


p = @(a,b) (N/t^3)*(beta0/pi)^(3/2)*exp(-beta0*(a.^2+b.^2)/t^2); %Boltzmann distribution

n = @(a,b) 1+p(a,b)*2*pi*alpha;
gradn= symfun(gradient(n,[a,b]),[a,b]);

ngradn =@(a,b) double( n(a,b).*gradn(a,b)); % represents n(r)*gradient(n(r)) is equal to acceleration function.



fcontour(@(x,y) log(n(x,y)),[-200,200,-200,200])%Plots contour map of n(r)

hold on
colorbar
% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 20;       
Tx0 = 1;      
Ty0 = 0; 
pos = [x0; y0; Tx0; Ty0];


w = trace2d(pos,1e-12,ngradn);

