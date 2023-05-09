%% Steps to computing the Navier-Stokes equations
% Step 1 1D linear convection


clc
clear all
close all

% Defining inital value of the variables

nx=101;                            % number of steps  in x
nt=50;                            % number of time steps    
dt=0.1;                           % step in time delta t
dx=10/(nx-1)                      % Step size
c=1;                              % Velocity   

% Establishing initial conditions that are known as 
% domain [0 10]
% u=2 @ 1<=x<=2
% u=1 @ everywhere else in [0 10]
% Boundary conditions u=1 @ x= 0 and x=10
% Discretizing using Forward Difference scheme for the time derivative 
% and the Backward Difference scheme for the space derivative.
u=ones(1,nx);
u(1,1/dx : 2/dx+1)=2
for p=1:1:nt
    un=u
    for i=2:1:nx
    
        u(i) = un(i) - c*dt/dx*(un(i)-un(i-1));
    end
    plot(linspace(0,10,nx),u)
    title('1D linear convection')
    xlabel('x');
    ylabel('u');
    pause(0.1)
end

%% Step 2 1D Convection

% Same code but change the constat velocity "c" to u(i)
% domain [0 10]
clc
clear all
close all


nx=101;                             % number of nodes in x
nt=50;                            % number of time steps    
dt=0.01;                           % step in time delta t
dx=10/(nx-1)                       % Step size

u=ones(1,nx);
u(1,1/dx : 2/dx+1)=2
%%% same discretization procedure as step 1
for p=1:1:nt
    un=u
    for i=2:1:nx
    
        u(i) = un(i) - u(i)*dt/dx*(un(i)-un(i-1));
    end
    plot(linspace(0,10,nx),u)
    title('1D non-linear Convection')
    xlabel('x')
    ylabel('u(i)')
    pause(0.1)
end

%% Step 3 1D Diffusion
% 
% domain [0 10]
clc
clear all
close all


nx=61;                              % number of nodes in x
nt=50;                              % number of time steps                  
dx=10/(nx-1);                       % Step size
vis=0.3;                            % viscosity
sigma=0.2;                          % Courant number
dt= sigma * dx^2/vis;               % step in time delta t
u=ones(1,nx);
u(1,1/dx : 2/dx+1)=2;
% Discretize the second-order derivative with a Central Difference a 
% combination of Forward Difference and Backward Difference of the first derivative.
for p=1:1:nt
    un=u
    for i=2:1:nx-1
   
        u(i) = un(i) + vis*dt/dx^2 * (un(i+1)-2*un(i)+un(i-1));
    end
    plot(linspace(0,10,nx),u)
    title('1D Diffusion')
    xlabel('x')
    ylabel('u')
    pause(0.1)
end

%% Step 4 Burgers equation in 1D

clc
clear all
close all
    
% Initial conditions and Boundary conditions are different than all steps
% ICs u=-2*dif*   and phi = exp(-(x-4*t)^2/(4*dif*(t+1))) + exp(-(x-4*t-2*pi)^2/(4*dif*(t+1)))
% BCs u(0)=u(2*pi)

syms x vis t
phi = exp(-(x-4*t)^2/(4*vis*(t+1))) + exp(-(x-4*t-2*pi)^2/(4*vis*(t+1)));
phiprime=diff(phi,x);
usym=-2*vis/phi*phiprime+4
ufunc=matlabFunction(usym)


% Time and local discretization

nx=101;                                % number of nodes in x
nt=100;                                % number of time steps 
vis=0.07;                             % viscosity
dx=2*pi / (nx-1);                     % Step size
dt=dx*vis;                            % step in time delta t
%timetrack = nt * dt;                  % time track for the graph
x=linspace(0,2*pi,nx);
t=0.1;
u=ufunc(vis,t,x);
un=zeros(1,nx);

% caclulation
% discretize sing forward difference for time, backward difference for 
% space and 2nd-order method for the second derivatives 
for n=1:1:nt
    un=u;
    for i=2:1:nx-1
        u(i) = un(i) - un(i) * dt/dx * (un(i)-un(i-1))+ vis * dt/dx^2 *(un(i+1)-2 *un(i) +un(i-1));
    end
    u(1) = un(1) - un(1) * dt/dx * (un(1)-un(nx-1))+ vis * dt/dx^2 *(un(2)-2 *un(1) +un(nx-1));
    u(nx)=u(1);
    plot(x,u)
   current_time = n * dt;
   title({'1D Burgers equation';['Time(\itt) = ', num2str(current_time), ' s']});
   xlabel('x')
   ylabel('y')
    pause(0.01)
    hold on
end
% this is the analytical solution which is in orange in the plot graph
u_analytical = ufunc(vis, n*dt,x);
plot(x, u_analytical,LineStyle="--",LineWidth=3);

%% Step 5 2D linear convection


clc
clear all
close all

% Defining inital value of the variables
% Domain in x [0 10] in y [0 10]
%setting inital parameters
nx=101;                            % number of nodes in x
ny=101;                            % number of nodes in y    
nt=100;                            % number of time steps
dx=10/(nx-1);                      % Step size in x
dy=10/(ny-1);                      % Step size in y
c=1;                               % Velocity  
sigma=0.2;                         % Courant number
dt=sigma * dx;                     % step in time delta t
timetrack = nt * dt;               % time track for the graph

%%%Definin the variables

x=linspace(0,10,nx);
y=linspace(0,10,ny);
u=ones(ny,nx);
un=ones(ny,nx);

u(2.5/dy:5/dy , 2.5/dx:5/dx)=2;

%meshing
[X,Y]=meshgrid(x,y);

% the timestep will be discretized as a forward difference and both spatial
% steps will be discretized as backward differences
for p=1:nt
    un=u
    for j=2:nx-1
        for i=2:ny-1
            u(i,j) = un(i,j) - c*dt/dx*(un(i,j)-un(i-1,j))-c*dt/dy*(un(i,j)-un(i,j-1));
            u(1:ny,1)=1;
            u(1,1:nx)=1;
            u(ny,1:nx);
            u(1:nx,ny)=1;
        end
    end
  
    surf(x,y,u);
    current_time = p * dt;
    title({'2D linear Convection';['Time(\itt) = ', num2str(current_time), ' s']});
    xlabel('(x)');
    ylabel('(y)');
    zlabel('u(x,y)')
    pause(0.01);
end

%% Step 6 2D Convection

clc
clear all
close all

% Defining inital value of the variables
% Domain in x [0 10] in y [0 10]
nx=101;                            % number of nodes in x
ny=101;                            % number of nodes in y    
nt=100;                            % number of time steps
dx=10/(nx-1);                      % Step size in x
dy=10/(ny-1);                      % Step size in y
sigma=0.2;
dt=sigma * dx;                     % step in time delta t
timetrack = nt * dt;               % time track

% Defining the variables

x=linspace(0,10,nx);
y=linspace(0,10,ny);
u=ones(ny,nx);
un=ones(ny,nx);
v=ones(ny,nx);
vn=ones(ny,nx);
u(2.5/dy:5/dy , 2.5/dx:5/dx)=2;

%meshing
[X,Y]=meshgrid(x,y);

%Finite difference calculation as before
for p=1:nt
    un=u;
    vn=v;
    for j=2:nx-1
        for i=2:ny-1
            u(i,j) = un(i,j) - un(i,j)*dt/dx*(un(i,j)-un(i-1,j))-vn(i,j)*dt/dy*(un(i,j)-un(i,j-1));
            v(i,j) = vn(i,j) - un(i,j)*dt/dx*(vn(i,j)-vn(i-1,j))-vn(i,j)*dt/dy*(vn(i,j)-vn(i,j-1));

            %%%boundary conditions
            u(1:ny,1)=1;
            u(1,1:nx)=1;
            u(ny,1:nx);
            u(1:nx,ny)=1;
            v(1:ny,1)=1;
            v(1,1:nx)=1;
            v(ny,1:nx);
            v(1:nx,ny)=1;
        end
    end
  
    surf(x,y,u);
    current_time = p * dt;
    title({'2D non-linear Convection';['Time(\itt) = ', num2str(current_time), ' s']});
    xlabel('(x)');
    ylabel('(y)');
    zlabel('u(x,y)')
    pause(0.01);
end

%% Step 7 2D diffusion
clc
clear all
close all

nx = 31;                                        % number of nodes in x
ny = 31;                                        % number of nodes in y
nt = 100;                                        % number of time steps
vis = 0.05;                                     % viscosity
dx = 2 / (nx - 1);                              % Step size in x
dy = 2 / (ny - 1);                              % Step size in y
sigma = 0.25;
dt = sigma * dx * dy / vis;                     % step in time delta t
timetrack = nt * dt;                            %time track for the graph

x=linspace (0,2,nx);
y=linspace(0,2,ny);

u=ones(ny,nx);
un=ones(ny,nx);

%ICs
u(0.5/dy:1/dy+1, 0.5/dx:1/dx+1)=2;

[X,Y]=meshgrid(x,y);
% discretizing forward difference in time and two second-order derivatives.
for p=1:nt+1
    un=u;
    for i=2:(nx-1)
        for j=2:(ny-1)
            %scheme of central differentiation
            u(i,j) = un(i,j) + vis * dt/dx^2 * (un(i+1,j) - 2*un(i,j) + un(i-1,j))+ vis * dt/dy^2 * (un(i,j+1)-2*un(i,j)+un(i,j-1));
        end
    end
    surf(X,Y,u);
    current_time = p * dt;
     title({'2D Diffusion';['Time(\itt) = ', num2str(current_time), ' s']});
    xlabel('(x)');
    ylabel('(y)');
    zlabel('Concentration')
    pause(0.05);
end

%% Step 8 2D Burgers equation

clc
clear all
close all

nx = 41;                                        % number of nodes in x
ny = 41;                                        % number of nodes in y
nt = 1000;                                       % number of time steps
vis = 0.05;                                     % viscosity coefficient
c=1                                             % velocity
dx = 2 / (nx - 1);                              % Step size in x
dy = 2 / (ny - 1);                              % Step size in y
sigma = 0.0009;
dt = sigma * dx * dy / vis;                     % step in time delta t
timetrack = nt * dt;                            %time track for the graph

x=linspace (0,2,nx);
y=linspace(0,2,ny);
u=ones(ny,nx);
v=ones(ny,nx);
un=u;
vn=v;
comb = ones(ny,nx)

u(0.5/dy:1/dy+1, 0.5/dx:1/dx+1)=2;
v(0.5/dy:1/dy+1, 0.5/dx:1/dx+1)=2;

[X,Y]=meshgrid(x,y);
for p=1:nt+1
    un=u;
    vn=v;
    for i=2:(nx-1)
        for j=2:(ny-1)
        u(i,j) = u(i,j) - (dt/dx) * u(i,j)*(u(i,j)-u(i-1,j)) - (dt/dy) * v(i,j)*(u(i,j)-u(i,j-1)) + (vis*dt/dx^2) * (u(i-1,j)-2*u(i,j)+u(i-1,j)) + (vis*dt/dy^2) * (u(i,j+1)-2*u(i,j)+u(i,j-1));    
        v(i,j) = v(i,j) - (dt/dx) * u(i,j)*(v(i,j)-v(i-1,j)) - (dt/dy) * v(i,j)*(v(i,j)-v(i,j-1)) + (vis*dt/dx^2) * (v(i-1,j)-2*v(i,j)+v(i-1,j)) + (vis*dt/dy^2) * (v(i,j+1)-2*v(i,j)+v(i,j-1)) ;
        u(1:ny,1)=1;
        u(1,1:nx)=1;
        u(ny,1:nx);
        u(1:nx,ny)=1;
        v(1:ny,1)=1;
        v(1,1:nx)=1;
        v(ny,1:nx);
        v(1:nx,ny)=1;
        end
    end
        surf(X,Y,u);
        pause(0.01);
end

%% Step 9 Laplace Equation

clc
clear all
close all


syms x y; % Define symbolic variables x and y

nx = 30;                                        % number of nodes in x
ny = 30;                                        % number of nodes in y
nt = 250;                                       % number of time steps
n_it=1000;                                      % number of iteration   
c=1                                             % velocity
dx = 2 / (nx - 1);                              % Step size in x
dy = 2 / (ny - 1);                              % Step size in y

x = 0:dx:2;                                     % domain
y = 0:dx:2;
p=zeros(nx,ny);
p(nx,:)=y

%discretized with central differences
for iit=1:n_it
    pn=p;
    for i =2:nx-1
        for j=2:ny-1
         p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2 + (pn(i,j-1)+pn(i,j+1))*dx^2)/(dx^2+dy^2)/2;
        end
    end
    p(2:nx-1,1) = p(2:nx-1,2);
    p(2:nx-1,ny) = p(2:nx-1,ny-1);
end

surf(x,y,p);
title('Numerical Solution of Laplace equation');

%%%% Analytical solution and the error
p_exact = zeros(size(x));
for n = 1:100
    p_exact = p_exact - 4/(n*pi)^2/sinh(2*pi*n)*sinh(n*pi*x).*cos(n*pi*y);
end
p_exact = x/4 + p_exact;

error=abs(p_exact-p);
disp(p_exact);
colormap(jet(256));
figure;
surf(x, y, error);
xlabel('x');
ylabel('y');
zlabel('Error');
title('Error between analytical and numerical solutions');


%% Step 10 Poisson's equation

clc
clear all
close all

nx = 30;                                         % number of nodes in x
ny = 30;                                         % number of nodes in y
nt = 100;                                        % number of time steps
xmin=0;
xmax=2;
ymin=0;
ymax=1;
n_it=1500;                                       % number of iteration   
dx = (xmax-xmin) / (nx - 1);                     % Step size in x
dy = (ymax-ymin) / (nx - 1);                     % Step size in y
p=zeros(nx,ny);
pn=zeros(nx,ny);
b=zeros(nx,ny);
x=xmin:dx:xmax;
y=ymin:dy:ymax;


b(floor(nx/4),floor(ny/4))=100;                  % source spike
b(floor(3*nx/4),floor(3*ny/4))=-100;             % source spike
[y,x]=meshgrid(y,x);
for nt=1:nt
    pn=p;
    for i=2:nx-1
        for j=2:ny-1
        p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2 + (pn(i,j-1)+pn(i,j+1))*dx^2 - b(i,j)*dx^2*dy^2)/(dx^2+dy^2)/2;   
        end
    end
p(2:nx-1,1) = p(2:nx-1,2);
p(2:nx-1,ny) = p(2:nx-1,ny-1);
surf(x,y,p)
colormap(jet(256))
pause(0.01)
end

%% Step 11 Final steps to solve the Navier-Stokes Equations

%initial conditions u,v,p=0 everywhere
%Boundar condition u=1 @ @y=2 ,   u,v=0 @ x=0,2 & y=0
% p=0 @ y=2  & dp\dx=0 @ x=0,2 & dp/dy=0 @y=0
clear all
close all

nx = 41;                                        % number of nodes in x
ny = 41;                                        % number of nodes in y
nt = 2000;                                        % number of time steps
xmin=0;
xmax=2;
ymin=0;
ymax=2
n_it=50;                                       % number of iteration   
dx = (xmax-xmin) / (nx - 1);                     % Step size in x
dy = (ymax-ymin) / (nx - 1);                            % Step size in y

x=linspace (0,2,nx);
y=linspace(0,2,ny);

[Y,X]=meshgrid(y,x);

rho=1;
vis=0.1;
dt=0.001;

%Preallocating
u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);
b=zeros(ny,nx);

%%%pressure field
for it=1:nt+1
    for i=2:(nx-1)
        for j=2:(ny-1)
            b(i,j) = rho * ((u(i+1,j) - u(i-1,j))/2/dx + (v(i,j+1) - v(i,j-1))/2/dy)/dt+((u(i+1,j) - u(i-1,j))/2/dx)^2 + 2*(u(i,j+1) - u(i,j-1))/2/dy*(v(i+1,j) - v(i,j-1))/2/dx + ((v(i,j+1) - v(i,j-1))/2/dy)^2;
        end
    end
    %%%poissons equation
for iit =1:n_it+1
    pn=p;
    for i=2:(nx-1)
        for j=2:(ny-1)
        p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2 + (pn(i,j+1)+pn(i,j-1))*dx^2 - b(i,j) * dx^2*dy^2)/(dx^2+dy^2)/2;
        end
    end
    p(1,:) = p(2,:);                           
    p(nx,:) = p(nx-1,:);                        
    p(:,1) = p(:,2);                            
    p(:,ny) = p(:,ny-1);                        
end

%%%%       velocity field     
un=u;
vn=v;
for i=2:nx-1
    for j=2:ny-1
        u(i,j) =un(i,j) - un(i,j) * dt / dx * (un(i,j) - un(i-1,j)) - vn(i,j) * dt / dy * (un(i,j) - un(i,j-1)) - 1/rho * (p(i+1,j) - p(i-1,j))*dt/2/dx + vis*dt/dx^2 * (un(i+1,j) - 2 * un(i,j) + un(i-1,j)) + vis*dt/dy^2 * (un(i,j+1) - 2 * un(i,j) + un(i,j-1));
        v(i,j) =vn(i,j) - vn(i,j) * dt / dy * (vn(i,j) - vn(i-1,j)) - vn(i,j) * dt / dy * (vn(i,j) - vn(i,j-1)) - 1/rho * (p(i,j+1) - p(i,j-1))*dt/2/dy + vis*dt/dy^2 * (vn(i+1,j) - 2 * vn(i,j) + vn(i-1,j)) + vis*dt/dy^2 * (vn(i,j+1) - 2 * vn(i,j) + vn(i,j-1));
    end
end
        u(1,:)  = 0;
        u(:,1)  = 0;
        u(nx,:) = 0;
        u(:,ny) = 1;   
        v(1,:)  = 0;
        v(:,ny) = 0;
        v(:,1)  = 0;
        v(nx,:) = 0;
        skip = 2;

end

figure;
contourf(X, Y, p, 'LineStyle', 'none');
colorbar;
    ylabel(colorbar, 'Pressure (p)');
    colormap(jet(256));
hold on;
contour(X, Y, p, 'LineWidth', 2, 'LineColor', 'non');
quiver(X(1:skip:end,1:skip:end), Y(1:skip:end,1:skip:end), u(1:skip:end,1:skip:end), v(1:skip:end,1:skip:end), 'Color', 'b', 'LineWidth', 1.5,MaxHeadSize=100);
hold off;
title('Fluid Flow Field in a Two-Dimensional Cavity');
xlabel('X');
ylabel('Y');
axis equal tight;


%% %% Step 12 Final simulation using Navier stokes equations


%%%% water simulation through a 2D square channel %%%%


%%% Initical conditions u,v,p= 0 everywhere
%%% boundary contitionsu u,v,p periodic @ x=0,2
%%%                               u,v=0 @ y=0,2
%%%                             dp/dy=0 @ y=0,2
clc
clear all
close all

nx = 21;                                        % number of nodes in x
ny = 21;                                        % number of nodes in y
nt = 50;                                        % number of time steps
xmin=0;
xmax=2;
ymin=0;
ymax=2
n_it=50;                                        % number of iteration   
dx = (xmax-xmin) / (nx - 1);                    % Step size in x
dy = (ymax-ymin) / (nx - 1);                    % Step size in y
F=1;
x=linspace (0,2,nx);
y=linspace(0,2,ny);

[Y,X]=meshgrid(y,x);

rho=1000;
vis=0.1;
dt=0.01;

u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);
b=zeros(ny,nx);

%pressure field

udiff=1;
stepcount=0;


while udiff >0.001 
    
    for i=2:(nx-1)
        for j=2:(ny-1)
            b(i,j) = rho * ((u(i+1,j) - u(i-1,j))/2/dx + (v(i,j+1) - v(i,j-1))/2/dy)+((u(i+1,j) - u(i-1,j))/2/dx)^2 - 2*(u(i,j+1) - u(i,j-1))/2/dy*(v(i+1,j) - v(i,j-1))/2/dx + ((v(i,j+1) - v(i,j-1))/2/dy)^2;
        end
    end
    %%% periodic term for x=0 what happens at x=0 happens at x=2
    for j=2:(ny-1)
        b(1,j) = rho * ((u(2,j) - u(nx,j))/2/dx + (v(i,j+1) - v(i,j-1))/2/dy)+((u(2,j) - u(nx,j))/2/dx)^2 - 2*(u(i,j+1) - u(i,j-1))/2/dy*(v(2,j) - v(1,j-1))/2/dx + ((v(i,j+1) - v(i,j-1))/2/dy)^2;
        end
    %%% periodic term for x=2 what happens at x=0 happens at x=2
    for j=2:(ny-1)
        b(nx,j) = rho * ((u(1,j) - u(nx-1,j))/2/dx + (v(i,j+1) - v(i,j-1))/2/dy)+((u(1,j) - u(nx-1,j))/2/dx)^2 - 2*(u(i,j+1) - u(i,j-1))/2/dy*(v(nx-1,j) - v(nx,j-1))/2/dx + ((v(i,j+1) - v(i,j-1))/2/dy)^2;
    
        %poissons equation
for iit =1:n_it+1
    pn=p;
    for i=2:(nx-1)
        for j=2:(ny-1)
        p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2 + (pn(i,j+1)+pn(i,j-1))*dx^2) /(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
        end
        %p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2 + (pn(i,j+1)+pn(i,j-1))*dx^2 - b(i,j) *dx^2 *dy^2)/(dx^2+dy^2)/2;
    end
    %%% periodic term for x=0
    for j=2:(ny-1)
         p(1,j)=((pn(2,j)+pn(nx,j))*dy^2 + (pn(1,j+1)+pn(1,j-1))*dx^2) /(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
    end
     %%% periodic term for x=2
    for j=2:(ny-1)
        p(nx,j)=((pn(2,j)+pn(nx-1,j))*dy^2 + (pn(nx,j+1)+pn(nx,j-1))*dx^2) /(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
        end
    p(1,:) = p(2,:);                           
    p(nx,:) = p(nx-1,:);                        
    p(:,1) = p(:,2);                                %dp/dy=0 at y=0
    p(:,ny) = p(:,ny-1);                            %dp/dy=0 at y=2 
end

%%%%       velocity field      %%
un=u;
vn=v;
for i=2:nx-1
    for j=2:ny-1
        u(i,j) =un(i,j) - un(i,j) * dt/dx * (un(i,j) - un(i-1,j)) - vn(i,j) * dt/dy * (un(i,j) - un(i,j-1))-dt / (2*rho*dx) * (p(i+1,j) - p(i-1,j)) + vis*(dt/dx^2 * (un(i+1,j) - 2 * un(i,j) + un(i-1,j)) + dt/dy^2 * (un(i,j+1) - 2 * un(i,j) + un(i,j-1))+ dt*F);
        v(i,j) =vn(i,j) - vn(i,j) * dt / dy * (vn(i,j) - vn(i-1,j)) - vn(i,j) * dt / dy * (vn(i,j) - vn(i,j-1))-dt / (2*rho*dx) * (p(i,j+1) - p(i,j-1)) + vis*(dt/dy^2 * (vn(i+1,j) - 2 * vn(i,j) + vn(i-1,j)) + dt/dy^2 * (vn(i,j+1) - 2 * vn(i,j) + vn(i,j-1)));
    end
end
%%%periodic for xx=0
for j=2:ny-1
        u(1,j) =un(1,j) - un(1,j) * dt / dx * (un(1,j) - un(nx,j)) - vn(1,j) * dt / dy * (un(1,j) - un(1,j-1))-dt / (2*rho*dx) * (p(2,j) - p(nx,j)) + vis*(dt/dx^2 * (un(2,j) - 2 * un(1,j) + un(nx,j)) + dt/dy^2 * (un(1,j+1) - 2 * un(1,j) + un(1,j-1))+ dt * F);
        v(1,j) =vn(1,j) - vn(1,j) * dt / dy * (vn(1,j) - vn(nx,j)) - vn(1,j) * dt / dy * (vn(1,j) - vn(1,j-1))-dt / (2*rho*dx) * (p(1,j+1) - p(1,j-1)) + vis*(dt/dy^2 * (vn(2,j) - 2 * vn(1,j) + vn(nx,j)) + dt/dy^2 * (vn(1,j+1) - 2 * vn(1,j) + vn(1,j-1)));
    end
%%%periodic for xx=0
for j=2:ny-1
        u(nx,j) =un(nx,j) - un(nx,j) * dt / dx * (un(nx,j) - un(nx-1,j)) - vn(nx,j) * dt / dy * (un(nx,j) - un(nx,j-1))-dt / (2*rho*dx) * (p(2,j) - p(nx,j)) + vis*(dt/dx^2 * (un(2,j) - 2 * un(nx,j) + un(nx-1,j)) + dt/dy^2 * (un(nx,j+1) - 2 * un(nx,j) + un(nx,j-1))+ dt * F);
        v(nx,j) =vn(nx,j) - vn(nx,j) * dt / dy * (vn(nx,j) - vn(nx-1,j)) - vn(nx,j) * dt / dy * (vn(nx,j) - vn(nx,j-1))-dt / (2*rho*dx) * (p(nx,j+1) - p(nx,j-1)) + vis*(dt/dy^2 * (vn(2,j) - 2 * vn(nx,j) + vn(nx-1,j)) + dt/dy^2 * (vn(nx,j+1) - 2 * vn(nx,j) + vn(nx,j-1)));
    end


        u(:,ny) = 0; 
        u(:,1)  = 0;
        v(:,1)  = 0;
        v(:,ny) = 0;
        
        
        
        
        udiff=(sum(sum(u))-sum(sum(un)))/sum(sum(u))
        stepcount=stepcount+1;
        quiver(x,y,u.',v.',1)
        colorbar;
        title('Flow through a 2D square channel ', ['stepcount ', num2str(stepcount)]);
        xlabel('X(cm)');
        ylabel('Y(cm)');
        axis equal tight
        hold off;
        pause(0.001)

    end
end



