%% FISHEYE LENS
%Sets up symbolic function environment
global ngradn
syms a b

%Input pressure function and symbolic computes refractive index and other
%functions

n = @(a,b) 2/(1+(a^2+b^2)); %Newton's fisheye lens refractive index
gradn= symfun(gradient(n,[a,b]),[a,b]);

ngradn =@(a,b) double( n(a,b).*gradn(a,b)); % represents n(r)*gradient(n(r)) is equal to acceleration function.


%Plots contour map of n(r)
fcontour(@(x,y) log(n(x,y)),[-2,2,-2,2])
hold on
colorbar


% trace() then computes the path of the ray.
%Initial conditions
x0 = -1/sqrt(2);       
y0 = 1/sqrt(2);       
Tx0 = 1;      
Ty0 = 0; 

w = trace(x0,y0,Tx0,Ty0);
figure
plot(w(:,4))

function output = trace(x0,y0,Tx0,Ty0)
tSpan = linspace(0,20,100);

      


% ode45 then solves this using the derivatives function. @derivatives calls
% the function derivatives as an object. This allows us just to change
% derivatives.
options = odeset('RelTol',1e-12,'Stats','off');
[t,w] = ode45(@ray_equation, tSpan, [x0; y0; Tx0; Ty0],options); %increase precision relative tolerance
    

% Visualising results. Possible animations and heat graphs?

%contourf(w,t)
plot(w(:,1),w(:,2),'-b');
hold on
plot(w(1,1),w(1,2),'bo','MarkerFaceColor','r');  % Mark the initial position with a red dot

axis equal
xlim([-5 5]);
ylim([-5 5]);

title('Path');
ylabel('y (m)');
xlabel('x (m)');


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