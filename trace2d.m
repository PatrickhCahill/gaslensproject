%% Plan to clean up trace
% Simplify trace to have optional inputs and save in external file. Inputs include
% 
% 1. Ray  position
% 2. Error tolerance
% 3. ngradn function
% 4. Create trace2d and trace 3d for different setups

%% Function

function output = trace2d(pos,tol,ngradn) %pos = [x0; y0; Tx0; Ty0]
tSpan = linspace(0,500,100);

      


% ode45 then solves this using the derivatives function. @derivatives calls
% the function derivatives as an object. This allows us just to change
% derivatives.
options = odeset('RelTol',tol,'Stats','off');
[t,w] = ode45(@ray_equation, tSpan, pos,options); %increase precision relative tolerance
    

% Visualising results. Possible animations and heat graphs?

%contourf(w,t)
plot(w(:,1),w(:,2),'-b');
hold on
plot(w(1,1),w(1,2),'bo','MarkerFaceColor','r');  % Mark the initial position with a red dot

axis equal
xlim([-200 200]);
ylim([-200 200]);

title('Path');
ylabel('y (m)');
xlabel('x (m)');


    % Ray equation for geometrical optics.
    function output = ray_equation(tf,wf)
        rf = [wf(1); wf(2)]; %radial vector as function of t
        Tf = [wf(3);wf(4)]; % optical ray vector as function of r.
        
        drdt = Tf;
        dTdt = ngradn(rf(1),rf(2));

        output = [drdt;dTdt];      

    end
    output = w;
end