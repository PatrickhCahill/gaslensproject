%% FISHEYE LENS CONVERGENCE TEST
%Sets up symbolic function environment
global ngradn
syms a b

%Input pressure function and symbolic computes refractive index and other
%functions

n = @(a,b) 2/(1+(a^2+b^2)); %Newton's fisheye lens refractive index
gradn= symfun(gradient(n,[a,b]),[a,b]);

ngradn =@(a,b) double( n(a,b).*gradn(a,b)); % represents n(r)*gradient(n(r)) is equal to acceleration function.

% trace() then computes the path of the ray.
%Initial conditions

errorlist = [];
for i = 1:1:12
    w = trace(i);
    error = abs(w(end,1).^2+w(end,2).^2-1);
    errorlist = [errorlist;10^(-1*i),error];
end
loglog(errorlist(:,1),errorlist(:,2))
set(gca, 'xdir','reverse');


function output = trace(tolerance)
tSpan = linspace(0,20,100);
x0 = -1;       
y0 = 0;       
Tx0 = 0;      
Ty0 = 1; 
      
tol = 10^(-1*tolerance);

% ode45 then solves this using the derivatives function. @derivatives calls
% the function derivatives as an object. This allows us just to change
% derivatives.
options = odeset('RelTol',tol,'Stats','off');
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