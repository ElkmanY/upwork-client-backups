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
N(1,3) = g;
A = [O, I; N, O];
B = [zeros(3,2); T];
C = eye(6);
D = zeros(6,2);

%% LQR gain
Q = 10000*diag([10 10 10 1 1 1]);
R = 1*diag([1 1]);
K = lqr(A,B,Q,R); % LQR gain
%% simulation parameters
h = 0.1;
tspan = 60;
simtime = 0:h:tspan;
nmax = length(simtime);
x0 = [0 500 0 0 -10 0]'; % initial condition
r = [0 0 0 0 0 0]'; % reference points
Umax = [0.5e3,44e3]';
Umin = [-0.5e3,-m*g]';
Xmax = [1000,1000, pi/2, 100, 100, pi/20]';
Xmin = [-1000, 0, -pi/2, -100, -100, -pi/20]';
%% simulation with linearized system
U = zeros(2,nmax);
X = zeros(6,nmax); 
X(:,1) = x0; 

for n = 2:nmax
	U(:,n-1) = max(min(10*K*(r-X(:,n-1)),Umax),Umin);
    d0 = zeros(3,1);
    d1 = 0.1*rand(3,1);
    X_dot = A*X(:,n-1) + B*U(:,n-1) + [zeros(3,1);d0];
    X(:,n) = max(min( X(:,n-1) + h*X_dot, Xmax),Xmin);
end
%%
figure
subplot(2,1,1)
plot(simtime,X(1,:),'k-');
hold on
plot(simtime,X(2,:),'k--');
hold off
legend({'x','y'})
ylabel('Position (m)');
subplot(2,1,2)
plot(simtime,X(3,:),'k-');
ylabel('Angle (rad)');
xlabel('Time (sec)');

figure
subplot(2,1,1)
plot(simtime,X(4,:),'k-');
hold on
plot(simtime,X(5,:),'k--');
hold off
legend({'x_{dot}','y_{dot}'})
ylabel('Velocity (m/sec)');
subplot(2,1,2)
plot(simtime,X(6,:),'k-');
ylabel('Angler velocity (rad/sec)');
xlabel('Time (sec)');
%%
figure
plot(simtime,U(1,:),'k-');
hold on
plot(simtime,U(2,:),'k--');
hold off
legend({'F_l','F_t'})
ylabel('Forces (N)');
xlabel('Time (sec)');