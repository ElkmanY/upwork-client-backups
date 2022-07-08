clear all
close all
%%
m = 20000; % mass - [kg]
J = 100000; %  [kg*m^2]
g = 1.6; % [m/s^2]
L = 4; % [m]
%% nonlinear system
M = diag([m,m,J]);
dx = @(x,u) [x(4:6); 
            u(1)*cos(x(3))/m + u(2)*sin(x(3))/m;
            -u(1)*sin(x(3))/m + u(2)*cos(x(3))/m;
            L/J*u(1)];
%% linearized system 
O = zeros(3,3); 
I = eye(3,3);
T = [1/m 0 L/J; 0 1/m 0]';
N = O;
N(1,3) = g/m;
A = [O, I; N, O];
B = [zeros(3,2); T];
%% LQR gain
Q = 1*diag([1 .01 1 1 1 1]);
R = 0.1*diag([.1 1]);
K = lqr(A,B,Q,R); % LQR gain
%% simulation parameters
h = 0.1;
tspan = 10000;
simtime = 0:h:tspan;
nmax = length(simtime);
% x0 = [-2 2 pi/2 0.1 0 0]'; % initial condition 1
x0 = [-20 20 0 0.1 2 0]'; % initial condition 2
r = [0 0 0 0 0 0]'; % reference points
Umax = [0.5e3,44e3]';
Umin = [-0.5e3,-m*g]';
dXmax = [1, .1, pi/6, 1, .1, pi/24]';
Xmax = [1000,1000, pi/2, 100, 100, pi/20]';
Xmin = [-1000, 0, -pi/2, -100, -100, -pi/20]';
%% simulation with linearized system
U = zeros(2,nmax);
X = zeros(6,nmax); 
X(:,1) = x0; 
for n = 2:nmax
	U(:,n-1) = max(min(K*(r-X(:,n-1)),Umax),Umin);
    X_dot = min(A*X(:,n-1) + B*U(:,n-1),dXmax);
    X(:,n) = max(min( X(:,n-1) + h*X_dot, Xmax),Xmin);
end
figure
plot(simtime,X(1:3,:));
%% simulation with nonlinear system
Un = zeros(2,nmax);
Xn = zeros(6,nmax); 
Xn(:,1) = x0; 
for n = 2:nmax
	Un(:,n-1) = max(min(K*(r-Xn(:,n-1)), Umax),Umin);
    X_dot = min(dx(Xn(:,n-1),Un(:,n-1)),dXmax);
    Xn(:,n) = max(min( Xn(:,n-1) + h*X_dot, Xmax),Xmin);
end
figure
plot(simtime,Xn(1:3,:));