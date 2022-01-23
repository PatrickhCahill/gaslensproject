%% BOLTZMANN LENS MULTIRAY
%Sets up symbolic function environment
global ngradn
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


p = @(a,b) (N/t^3)*(beta0/pi)^(3/2)*exp(-beta0*(a^2+b^2)/t^2); %Boltzmann distribution

n = @(a,b) 1+p(a,b)*2*pi*alpha;
gradn= symfun(gradient(n,[a,b]),[a,b]);

ngradn =@(a,b) double( n(a,b).*gradn(a,b)); % represents n(r)*gradient(n(r)) is equal to acceleration function.



% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 5;       
Tx0 = 1;      
Ty0 = 0; 
for i = 0:5:45
    w = trace(x0,y0+i,Tx0,Ty0);
    final = extrapolate(w);
    plot(final(:,1),final(:,2),'-b')
    hold on
end
xlim([0,20e9])
xlabel("x (m)")
ylabel("y (m)")

function output = extrapolate(w)
    x1 = w(99,1);
    x2 = w(100,1);
    y1 = w(99,2);
    y2 = w(100,2);
        
    grad = (y2-y1)/(x2-x1);
    xend = (-y2)/grad+x2;
    output = [x2,y2;xend,0];
end
function output = trace(x0,y0,Tx0,Ty0)
tSpan = linspace(0,500,100);

      


% ode45 then solves this using the derivatives function. @derivatives calls
% the function derivatives as an object. This allows us just to change
% derivatives.
options = odeset('RelTol',1e-12,'Stats','off');
[t,w] = ode45(@ray_equation, tSpan, [x0; y0; Tx0; Ty0],options); %increase precision relative tolerance
    
    % Ray equation for geometrical optics.
    function output = ray_equation(tf,wf)
        global ngradn
        rf = [wf(1); wf(2)]; %radial vector as function of t
        Tf = [wf(3);wf(4)]; % optical ray vector as function of r.
        
        drdt = Tf;
        dTdt = ngradn(rf(1),rf(2));

        output = [drdt;dTdt];      

    end
    output = w;
end