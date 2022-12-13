% Simulation of the generic 3-state system with internal model
% form for reference input tracking
% the following parameters are used in the model function file
% for simulation
clear all

global K A B C D ro

%SETTING VARIABLES
g = 9.81; %m/s^2
L = 11.2e-2; %m
mp = 381e-3; %kg
Ip = 0.00616; %kg m2
mw = 36e-3; %kg
i_cmw = 7.46e-6; %kg m^2
rw = 2.1e-2; %m
R = 4.4; %ohm
kb = 0.495; %Vs/rad
kt = 0.470; %Nm/A
q1 = L*mp;
q2 = Ip + (L^2)*mp;
q3 = (kb*kt)/R;
q4 = -mp - mw - i_cmw/rw^2;

%SETTING MATRICES
A12 = (g*q1*q4)/(q1^2+q2*q4);
A22 = (q3 * (q1- q4*rw))/((q1^2+q2*q4)*rw);
A42 = (q3*(-q1 + q4*rw))/((q1^2+q2*q4)*(rw^2));
A14 = (g*q1^2)/(q1^2+q2*q4);
A24 = -(q3 * (q2+ q1*rw))/((q1^2+q2*q4)*rw);
A44 = (q3*(q2 + q1*rw))/((q1^2+q2*q4)*(rw^2));
B2 = (kt*(q1-q4*rw))/(R*(q1^2+q2*q4)*rw);
B4 = 20.6;
% B4 = (kt*(-q2+q1*rw))/(R*(q1^2-q2*q4)*rw);

%SETTING UP STATE SPACE MATRIX
A = [0, 1, 0, 0; 
    A12, A22, 0 , A42;
    0, 0, 0, 1;
    A14, A24, 0 , A44];

B = [0;B2;0;B4];

% C = [1, 0, 0, 0; 0, 0, 1, 0]; %C matrix to isolate alpha and x
C = [1,0,1,0];
D = 0;

% This design uses an internal model (IM)
% for IM form (see eq. 11.81 in Dorf & Bishop)
Aim = [0,C;zeros(4,1),A]; Bim = [0;B(:,1)];
% NOTE that the IM only uses the 1st column of B for the control input
% The second input is a disturbance
% form internal model controllability matrix
Pimc = ctrb(Aim,Bim);
% confirm system is completely controllable; det(Pc) not zero
det(Pimc)

% Select poles for IM
p1 = -0.391; p2 = -3.691 + 1.968i; p3 = -3.691 - 1.968i; p4 = -1064; p5 = -1000;
% poles = [p1, p2, p3, p4, p5];

Q = [1,0,0,0;
        0,0.001,0,0;
        0,0,100,0;
        0,0,0,1]; 
R = 1;
[K_LQR,S_LQR,P_LQR] = lqr(A,B,Q,R);
% select gains
poles = [P_LQR',min(P_LQR)-10]
K = place(Aim,Bim,poles );


% Set up for simulation using RK4
% define the initial states of the plant
% make sure they are compatible
xo = [deg2rad(1);0;0;0]; % plant initial states
yo = C*xo; % initial output
ro = yo; % let the reference initialize equal to yo
eo = yo - ro; % initial error state
uo = 0; % initial control input

% initialize z; z = xdot, so zo comes from xdoto
zo = A*xo+B(:,1)*uo;
% construct initial condition vector for IM system; include plant
xco = [eo;zo;xo;uo];
% 2(n+1) = 8 states
% Simulation parameters
dt = 0.00002; % time interval for fixed-step simulation
t0 = 0.0; % start time
tf = 1; % final time
N = floor((tf-t0)/dt); % number of steps

% call RK4
[tc,xc] = rk4fixed(@four_state_IM,[t0 tf],xco,N);
% Compute defined outputs
for i = 1:N
[xd,y] = four_state_IM(tc(i),xc(i,:));
e(i) = y(1);
r(i) = y(2);
end

% Plotting
figure(2)
subplot(3,1,1), plot(tc,e), legend('e = y-r')
title('Error state'), grid on
subplot(3,1,2), plot(tc,xc(:,10)), legend('uc')
title('Control input'), grid on
subplot(3,1,3), plot(tc,r,tc,xc(:,7),tc,xc(:,9)),legend('r','\alpha','x')
title('r, \alpha, x'), grid on
% % ylim([0,1.5])
% subplot(4,1,4), plot(tc,xc(:,2)), legend('x 2')
% ylim([-0.5,0.5])
