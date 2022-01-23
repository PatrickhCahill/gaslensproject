%% 3D BOLTZMANN LENS - Hexagonal Distribution
%Sets up symbolic function environment
global ngradn
syms x y z

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
radius = 25;
beta0 = particlemass/(2*k*temp0);

%Input pressure function and symbolic computes refractive index and other
%functions


pcentred = @(x,y,z) (N/t^3)*(beta0/pi)^(3/2).*exp(-beta0.*(x.^2+y.^2+z.^2)./t^2); %Boltzmann distribution



p = @(x,y,z) pcentred(x,y-radius*cos(pi/3),z-radius*sin(pi/3))+pcentred(x,y-radius*cos(2*pi/3),z-radius*sin(2*pi/3))+pcentred(x,y-radius*cos(3*pi/3),z-radius*sin(3*pi/3))+pcentred(x,y-radius*cos(4*pi/3),z-radius*sin(4*pi/3))+pcentred(x,y-radius*cos(5*pi/3),z-radius*sin(5*pi/3))+pcentred(x,y-radius*cos(6*pi/3),z-radius*sin(6*pi/3));
n = @(x,y,z) 1+p(x,y,z)*2*pi*alpha;
gradn= symfun(gradient(n,[x,y,z]),[x,y,z]);

ngradn =@(x,y,z) double( n(x,y,z).*gradn(x,y,z)); % represents n(r)*gradient(n(r)) is equal to acceleration function.


% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 2;
z0 = 2;
Tx0 = 1;      
Ty0 = 0; 
Tz0 = 0;


w = trace(x0,y0,z0,Tx0,Ty0,Tz0);
ns = log(n(w(:,1)',w(:,2)',w(:,3)')-1);
surf([w(:,1)';w(:,1)'],[w(:,2)';w(:,2)'],[w(:,3)';w(:,3)'],[ns;ns],'facecolor','none','edgecolor','interp','linewidth',2);
cb = colorbar;
cb.Label.String = '{log (n(r)-1)}';
axis equal
xlim([-200 200]);
ylim([-200 200]);
zlim([-200 200]);

title('Path');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
hold on
for i = -2:2
    for j = -2:2
        if i ~= 0 || j ~=0
            x0 = -180;       
            y0 = 2+i*50;
            z0 = 2+j*50;
            Tx0 = 1;      
            Ty0 = 0; 
            Tz0 = 0;
            w = trace(x0,y0,z0,Tx0,Ty0,Tz0);
            ns = log(n(w(:,1)',w(:,2)',w(:,3)')-1);
            surf([w(:,1)';w(:,1)'],[w(:,2)';w(:,2)'],[w(:,3)';w(:,3)'],[ns;ns],'facecolor','none','edgecolor','interp','linewidth',2);
        end
    end
end




function output = extrapolate(w)
    x1 = w(99,1);
    x2 = w(100,1);
    y1 = w(99,2);
    y2 = w(100,2);
        
    grad = (y2-y1)/(x2-x1);
    xend = (-y2)/grad+x2;
    output = [x2,y2;xend,0];
end
function output = trace(x0,y0,z0,Tx0,Ty0,Tz0)
tSpan = linspace(0,500,100);

      


% ode45 then solves this using the derivatives function. @derivatives calls
% the function derivatives as an object. This allows us just to change
% derivatives.
options = odeset('RelTol',1e-12,'Stats','off');
[t,w] = ode45(@ray_equation, tSpan, [x0; y0; z0; Tx0; Ty0; Tz0],options); %increase precision relative tolerance
    

% Visualising results. Possible animations and heat graphs?

%contourf(w,t)


    % Ray equation for geometrical optics.
    function output = ray_equation(tf,wf)
        global ngradn
        rf = [wf(1); wf(2); wf(3)]; %radial vector as function of t
        Tf = [wf(4); wf(5); wf(6)]; % optical ray vector as function of r.
        
        drdt = Tf;
        dTdt = ngradn(rf(1),rf(2),rf(3));

        output = [drdt;dTdt];      

    end
    output = w;
end