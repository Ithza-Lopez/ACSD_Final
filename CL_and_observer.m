clear all

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

% %desired poles taken from LQR in previous code
% Q = [1,0,0,0;
%         0,0.001,0,0;
%         0,0,100,0;
%         0,0,0,1]; 
% R = 1;
% [K_LQR,S_LQR,P_LQR] = lqr(A,B,Q,R);
% K = place(A,B,P_LQR);
% p1 = -0.391; p2 = -3.691 + 1.968i; p3 = -3.691 - 1.968i; p4 = -1064;
% p1 = -0.391; p2 = -3.691 + 1.968i; p3 = -3.691 - 1.968i; p4 = -0.04;
% p1 = -1; p2 = -2; p3 =  -3; p4 = -4;
% poles = [p1, p2, p3, p4];
x0 = [deg2rad(20), 0, 0, 0];
e0 = [deg2rad(1), 0.0, 0.0, 0.0];
tlim = 0:0.001:1;
% check controllability
Pc=ctrb(A,B);
if rank(Pc)==length(A)
%     Q = eye(4)*100;
    Q = [1,0,0,0;
        0,0.001,0,0;
        0,0,100,0;
        0,0,0,1]; 
    R = 1;
    [K_LQR,S_LQR,P_LQR] = lqr(A,B,Q,R);
    K = place(A,B,P_LQR);
    % compare response of OL system and state-feedback regulator
    Ac = A-B*K;
    Bc = [0;0;0;0];
    sys = ss(A,B,C,D);
    sysc = ss(Ac,Bc,C,D);
    [y,t,x] = initial(sys, x0 ,tlim);
    [yc,tc,xc] = initial(sysc,x0,tlim);
    
    figure(1)
    subplot(2,1,1), plot(t,x(:,1),t,x(:,3)), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
    title('OL response');
    subplot(2,1,2), plot(tc,xc(:,1),tc,xc(:,3)), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
    title('CL response')
    hold on
    grid on
    fsctrb = 1;
    
    figure(2)
    subplot(2,1,1)
    plot(tc,xc(:,1),tc,xc(:,3)), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
    grid on
    title('Close-Loop Response')
    subplot(2,1,2)
    rlocus(sysc)
    grid on
    title('Root Locus CL')
    hold off
    
    
else
    sprintf('not controllable')
    fsctrb = 0;
end

% check observability
Po=obsv(A,C);
if rank(Po)==length(A)
    dom_pole = max(real(P_LQR));
    pole_muliplier = dom_pole*10;
    po = [pole_muliplier, pole_muliplier-1, pole_muliplier-2, pole_muliplier-3];
    L = place(A',C',po)';
    
    % build observer system
    Ao = A-L*C;

    sysoe = ss(Ao,[0;0;0;0],eye(4),[0;0;0;0]);
    [yoe,toe,xoe]=initial(sysoe,x0,tlim);
    figure(3)
    plot(toe,xoe(:,1),toe,xoe(:,3)), legend({'e $\alpha$', 'e $x$'}, 'Interpreter', 'latex')
    grid on
    fsobsv = 1;
else
    sprintf('not observable')
    fsobsv = 0;
end

if (fsctrb==1)&&(fsobsv==1)
%     build the integrated system
    Aa = [Ac,B*K;zeros(4),Ao];
    Ba = [B;0;0;0;0];
    Ca = [0,1,0,1,0,0,0,0];
    Da = 0;
    sysa = ss(Aa,Ba,Ca,Da);
    [ya,ta,xa]=initial(sysa,[x0,e0],tlim);
    
    figure(4)
    
    plot(tc,xc(:,1),'--r',tc,xc(:,3),'--b')
    hold on
    plot(ta,xa(:,1),'-r',ta,xa(:,3),'-b')
    title('Comparing CL and observer')
    legend({'CL $\alpha$', 'CL $x$', 'Obs $\alpha$', 'Obs $x$'}, 'Interpreter', 'latex')
    grid on
%     subplot(2,1,2), plot(ta,xc(:,1)-xa(:,1),ta,xc(:,3)-xa(:,3)), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
%     title('Difference in responses')
%     grid on
%     hold off
%     subplot(3,1,1), plot(tc,xc(:,1),':r',tc,xc(:,3),':b'), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
%     title('controlled with FSFB')
%     grid on
%     subplot(3,1,2), plot(ta,xa(:,1),'-r',ta,xa(:,3),'-b'), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
%     title('controlled with observer FSFB')
%     grid on
%     subplot(3,1,3), plot(ta,xc(:,1)-xa(:,1),ta,xc(:,3)-xa(:,3)), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
%     title('Difference in responses')
%     grid on
%     hold off
    
    
    figure(5)
    title('Voltage output')
    VoltK = -K*xc'; VoltObs = -K*xa(:,1:4)';
    plot(t, VoltK); hold on
    plot(ta, VoltObs);
    legend('Controlled Voltage', 'Voltage with Observer')
    hold off
    grid on
    
%     figure(6) %OP reponse and root locus
%     subplot(2,1,1)
%     plot(t,x(:,1),t,x(:,3)), legend({'$\alpha$', '$x$'}, 'Interpreter', 'latex')
%     grid on
%     title('Open-Loop Response')
%     subplot(2,1,2)
%     rlocus(sys)
%     grid on
%     title('Root Locus OP')
%     hold on
end

%reference input

