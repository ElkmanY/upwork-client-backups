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
Q = 20000*diag([10 10 10 1 1 1]);
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
dUmax = 0.1*Umax;
Xmax = [1000,1000, pi/6, 10, 10, pi/60]';
Xmin = [-1000, 0, -pi/6, -10, -10, -pi/60]';
%% simulation with linearized system
U = zeros(2,nmax);
X = zeros(6,nmax); 
Xm = zeros(6,nmax); 
X(:,1) = x0; 

for n = 2:nmax
    d = [0.1*rand(2,1);0];
    d = zeros(3,1);
    d_hat = m*d(1:2);
	U(:,n-1) = max(min(10*K*(r-Xm(:,n-1)) + d_hat,Umax),Umin);  
    if n>2
        dU = min(max(U(:,n-1) - U(:,n-2),-dUmax),dUmax);
        U(:,n-1) = U(:,n-2) + dU;
    end
    X_dot = A*X(:,n-1) + B*U(:,n-1) + [zeros(3,1);d];
    X(:,n) = max(min( X(:,n-1) + h*X_dot, Xmax),Xmin);
    
    x = X(1,n); y = X(2,n);
    Am = atan2(x+100,y) + h*h*(0.06*rand(1)-0.03);
    Bm = atan2(x-100,y) + h*h*(0.06*rand(1)-0.03);
    Xm(1,n) = 100*(tan(Am)+tan(Bm))/(tan(Am)-tan(Bm));
    Xm(2,n) = 200/(tan(Am)-tan(Bm));
    Xm(6,n) = X(6,n) + h*(1*rand(1)-0.5);
    Xm(3,n) = Xm(3,n-1) + (Xm(6,n)-Xm(6,n-1));
    Xm(4:5,n) = Xm(4:5,n-1) + h*(X_dot(4:5) + h*(0.2*rand(2,1)-0.1));
    
    Xm(:,n) = max(min( Xm(:,n), Xmax),Xmin);
end
%%
figure
subplot(2,1,1)
plot(simtime,Xm(1,:),'k-');
hold on
plot(simtime,Xm(2,:),'k--');
hold off
legend({'x','y'})
ylabel('Position (m)');
subplot(2,1,2)
plot(simtime,Xm(3,:),'k-');
ylabel('Angle (rad)');
xlabel('Time (sec)');

figure
subplot(2,1,1)
plot(simtime,Xm(4,:),'k-');
hold on
plot(simtime,Xm(5,:),'k--');
hold off
legend({'x_{dot}','y_{dot}'})
ylabel('Velocity (m/sec)');
subplot(2,1,2)
plot(simtime,Xm(6,:),'k-');
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
