    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Jalpesh Bhadra                  %
%            EEE 588 Project                 %
%       Control of two dof robotic arm       %

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters
syms L1 L2 q1 q2
% L1 = Length of Link 1
% L2 = Length of Link 2
% q1 = Angle of Link 1
% q2 = Angle of Link 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation of Forword Kinematics
%T01 = Transition matrix from Base frame to Frame 1
T01 = [cos(q1) -sin(q1) 0 L1*cos(q1);
        sin(q1) cos(q1) 0 L1*sin(q1);
        0            0  0     0     ;
        0             0  0     1     ];
    
%T12 = Transition matrix from frame 1 to frame 2
T12 = [cos(q2) -sin(q2) 0 L1*cos(q2);
        sin(q2) cos(q2) 0 L1*sin(q2);
        0            0  0     0     ;
        0             0  0     1     ];
    
 T02 = simplify(T01*T12);
 
%  T02 =
%  
% [ cos(q1 + q2), -sin(q1 + q2), 0, L1*(cos(q1 + q2) + cos(q1))]
% [ sin(q1 + q2),  cos(q1 + q2), 0, L1*(sin(q1 + q2) + sin(q1))]
% [            0,             0, 0,                           0]
% [            0,             0, 0,                           1]
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %      Inverse Kinematic 
%  syms xe ye 
%  % xe = X coordinate of end effector
%  % ye = Y coordinate of end effector
%  con = (xe^2 + ye^2 - L1^2 - L2^2)/(2*L1*L2)
%  cosq2 = con
%  sinq2 = (1-(cosq2)^2)^0.5
%  q2 = atan2(sinq2,cosq2)
%  q1 = atan2(ye/xe) - atan2((L2*sin(q2))/(L1+L2*cos(q2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Parameters  
syms a1 a2 L1 L2 IL1 IL2 mL1 mL2 kr1 kr2 mm1 mm2 Im1 Im2 T1 T2 x1 x2 x3 x4
% x1 = tetha 1 (Joint angle 1)
% x2 = tetha 2 (Joint angle 2)
% x3 = tetha 1 dot (Angular velocity of joint 1)
% x4 = tetha 2 dot (Angular velocity of joint 2)
a1=1;a2=1; %Length of Link 1 is a1 (m) ; Lenght of Link 2 is a2 (m)
L1=0.5; L2=0.5; %Distance of CG on link - L1 (m); L2(m)
mL1=50; mL2=50; %mL1 and mL2 are weight (in kg) of link 1 and 2 respectively 
IL1=10; IL2=10; %Il1 and Il2 is moment of inertia (kg.mm^2) of link 1 and 2 respectively
kr1=100; kr2=100; %kr1, kr2 is gear ratio of motor 1 and 2 respectively. 
mm1=5; mm2=5; %mm1 and mm2 is weight (in kg) of motor 1 and 2 respectively
Im1=0.01; Im2=0.01; %Im1 and Im2 is moment of inertia for motor 1 and 2 respectively
%T1 and T2 is torque applied to motor 1 and 2 respectively 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Dynamics of 2 link dof robot
%Mass/Inertia matrix M is a 2x2 matrix
M(1,1) = IL1+ mL1*L1^2 + kr1^2*Im1 + IL2 + mL2 *(a1^2+L2^2+2*a1*L2*cos(x2)) + Im2 + mm2*a1^2;
M(1,2) = IL2 + mL2*(L2^2+a1*L2*cos(x2))+ kr2*Im2;
M(2,1) = M(1,2);
M(2,2) = IL2 + mL2*L2^2 + kr2^2*Im2;

%Centrifugal/Corlois matrix C is a 2x2 matrix
syms h g x1 x2 x3 x4 
h = -mL2*a1*L2*sin(x2);
C(1,1) = h*x4;
C(1,2) = h*(x3+x4);
C(2,1) = -h*(x3);
C(2,2) = 0;

%Gravity matrix G
g = 9.81; %(m/s^2)
G(1,1) = (mL1*L1 + mm2*a1 + mL2*a1)*g*cos(x1) + mL2*L2*g*cos(x1+x2);
G(1,2) = mL2*L2*g*cos(x1+x2);
G(2,1) = G(1,2);
G(2,2) = G(1,2);

% Desription of Non linear model
 % M*q_double_dot + C(q,q_dot) + G(q) = T 
 % q = Joint angles 
 % q_dot = Joint velocities
 % q_double_dot = Joint acceleration
 % M = Mass/Inertia Matrix
 % C = Cemtrifugal/Corolis matrix \\
 % G = gravity matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearizing the model
% Linearizing about tetha 1 = pi/2 rad and tetha 2 =0 ;
% State vector x = [x1 x2 x3 x4]
% x1 = tetha 1 (Joint angle 1)
% x2 = tetha 2 (Joint angle 2)
% x3 = tetha 1 dot (Angular velocity of joint 1)
% x4 = tetha 2 dot (Angular velocity of joint 2)
% xequilibrium = x1_eq, x2_eq, x3_eq, x4_eq
x1_eq = pi/2;
x2_eq = 0;
x3_eq = 0;
x4_eq = 0;

%% Linearize mass/inertia matrix
Meq(1,1) = IL1+ mL1*L1^2 + kr1^2*Im1 + IL2 + mL2 *(a1^2+L2^2+2*a1*L2*cos(x2_eq)) + Im2 + mm2*a1^2;
Meq(1,2) = IL2 + mL2*(L2^2+a1*L2*cos(x2_eq))+ kr2*Im2;
Meq(2,1) = Meq(1,2);
Meq(2,2) = IL2 + mL2*L2^2 + kr2^2*Im2;

Geq(1,1) = -(mL1*L1 + mm2*a1 + mL2*a1)*g*sin(x1_eq) - mL2*L2*g*sin(x1_eq+x2_eq);
Geq(1,2) = -mL2*L2*g*sin(x1_eq+x2_eq);
Geq(2,1) = Geq(1,2);
Geq(2,2) = Geq(1,2);

Meqinverse = inv(Meq);
A=[zeros(2,2) eye(2,2) ; -Meqinverse*Geq zeros(2,2)];
B=[zeros(2,2) ; Meqinverse]; 
C = [eye(2,2) zeros(2,2)];
D = zeros(2,2);
% 
% A =
% 
%          0         0    1.0000         0
%          0         0         0    1.0000
%     4.0421    0.6419         0         0
%     0.4017    1.7479         0         0

% B =
% 
%          0         0
%          0         0
%     0.0043   -0.0017
%    -0.0017    0.0088
% C =
% 
%      1     0     0     0
%      0     1     0     0
% 
% 
% D =
% 
%      0     0
%      0     0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scaling 
r2d = 180/pi;
sx  = diag( [ r2d, r2d, r2d, r2d ] );   % Convert radians to degrees
sy  = diag( [ r2d, r2d] );              % Convert radians to degrees

A=sx*A*inv(sx);
B=sx*B;
C=sy*C*inv(sx);

% % A =
% % 
% %          0         0    1.0000         0
% %          0         0         0    1.0000
% %     4.0421    0.6419         0         0
% %     0.4017    1.7479         0         0
% % 
% % 
% % B =
% % 
% %          0         0
% %          0         0
% %     0.2482   -0.0983
% %    -0.0983    0.5066
% % 
% % 
% % C =
% % 
% %      1     0     0     0
% %      0     1     0     0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Analysis of system

% Number of states = ns
% Number of controls = nc
[ns,nc]=size(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Modal Analysis of System
%   Computing eigenvalues and eigenvectors for system
[V,Dia] = eig(A);  %V contains eigenvector matrix
                %D contains eigenvalue matrix
                
% % V =
% % 
% %    -0.4346    0.1589    0.4346   -0.1589
% %    -0.0727   -0.5945    0.0727    0.5945
% %     0.8854   -0.2035    0.8854   -0.2035
% %     0.1481    0.7615    0.1481    0.7615
% % 
% % 
% % Dia =
% % 
% %    -2.0370         0         0         0
% %          0   -1.2808         0         0
% %          0         0    2.0370         0
% %          0         0         0    1.2808

% The right half plane pole at s = 2.0370 is a standard inverted pendulum
% instability. It is associated with theta1_dot
% 
% The left half plane pole at s = -2.0370 is a standard inverted pendulum damping mode
% It is associated with theta1_dot
% 
% The right half plane pole at s = 1.2808 is a standard inverted pendulum
% instability. It is associated with theta2_dot
% 
% The left half plane pole at s = -1.2808 is a standard inverted pendulum damping mode
% It is associated with theta2_dot

%Let us compute natural modes (temdencies) of system

tinit = 0;
tinc = 0.001;
tfin = 0.3;
t=[tinit:tinc:tfin]; %Uniform time vector
u =[0*t' 0*t']; %Zero input response

% Exciting mode at s = 2.0370 which is a fast instability
% This mode is associated with a pole at s= 2.0370 and with theta1_dot
x = lsim(ss(A,B, eye(ns,ns), 0*ones(ns,nc)),u,t,V(:,3));
figure(1);plot(t,x);
grid
title('Fast Instability: x_o = [0.4346 0.0727 0.8854 0.1481 ]^T')
ylabel('States (deg, deg/sec)')
xlabel('Time (seconds)')


% Exciting mode at s= 1.2080 which is slow instability.
% This mode is associated with a pole at s=1.2080 and theta2_dot
x = lsim(ss(A,B, eye(ns,ns), 0*ones(ns,nc)),u,t,V(:,4));
figure(2);plot(t,x);
grid
title('Slow Instability: x_o = [-0.1589 0.5945 -0.2035 0.7615 ]^T')
ylabel('States (deg, deg/sec)')
xlabel('Time (seconds)')


% Exciting mode at s = -2.0370 which is a stable
% This mode is associated with a pole at s= -2.0370 and with theta1_dot
x = lsim(ss(A,B, eye(ns,ns), 0*ones(ns,nc)),u,t,V(:,1));
figure(3);plot(t,x);
grid
title('quick stable mode: x_o = [-0.4346 -0.0727 0.8854 0.1481]^T')
ylabel('States (deg, deg/sec)')
xlabel('Time (seconds)')

% Exciting mode at s= 1.2080 which is slow instability.
% This mode is associated with a pole at s=1.2080 and theta2_dot
x = lsim(ss(A,B, eye(ns,ns), 0*ones(ns,nc)),u,t,V(:,2));
figure(4);plot(t,x);
grid
title('slow stable mode: x_o = [0.1589 -0.5945 -0.2035 0.7615]^T')
ylabel('States (deg, deg/sec)')
xlabel('Time (seconds)')

%% Damping and natural frequency for each pole
sys = ss(A,B, eye(ns,ns), 0*ones(ns,nc),0.5);
damp(sys) %damp gives many output collect output in vectors as done in F8
%tau = 1./(wn.*zeta) 

% % Step response charcteristic
% S = stepinfo(sys)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%       Transmission Zeros
plantzeros = tzero(sys);
 % System has no finite transmission zeros
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       System Transfer Functions
%plat_tf = zpk (sys)
plant_tf = zpk(ss(A,B,C,D));
% 
% plat_tf =
%  
%   From input 1 to output...
%              0.24824 (s-1.415) (s+1.415)
%    1:  ---------------------------------------
%        (s+2.037) (s+1.281) (s-1.281) (s-2.037)
%  
%             -0.098283 (s-2.249) (s+2.249)
%    2:  ---------------------------------------
%        (s+2.037) (s+1.281) (s-1.281) (s-2.037)
%  
%   From input 2 to output...
%             -0.098283 (s-2.249) (s+2.249)
%    1:  ---------------------------------------
%        (s+2.037) (s+1.281) (s-1.281) (s-2.037)
%  
%               0.50663 (s-2.03) (s+2.03)
%    2:  ---------------------------------------
%        (s+2.037) (s+1.281) (s-1.281) (s-2.037)
%  
% Continuous-time zero/pole/gain model.
w=logspace(-1,2,500);
figure(5);
subplot(2,2,1)
bode(plant_tf(1,1),w)
title('Output \theta 1 Input \tau 1')

subplot(2,2,2)
bode(plant_tf(1,2),w)
title('Output \theta 2 Input \tau 1')

subplot(2,2,3)
bode(plant_tf(2,1),w)
title('Output \theta 1 Input \tau 2')

subplot(2,2,4)
bode(plant_tf(2,2),w)
title('Output \theta 2 Input \tau 2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%       Controllability 
%
contro = [B A*B A*A*B A*A*A*B]; %Controlability matrix
rcontro = rank(contro);
% system is controllable; rank os controlability matrix is 4
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%       Observability
observ = [C; C*A; C*A*A; C*A*A*A];
robserv = rank(observ);
% System is obervable; rank of observability matrix is 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% PLANT STEP RESPONSE
%
t = [0:0.1:600];
x = step(ss(A, B, eye(4,4), 0*ones(4,2)), t); % Step response for each state and each input

%
% Step Response due to \tau 1
%
%\theta 1 response
figure(6)
plot(t,x(:,1,1))
grid
title('\theta 1 response for \tau 1 = Unit Step')
ylabel('\theta 1 (deg)')
xlabel('Time (seconds)')
%return
%
%\theta 2 response
figure(7)
plot(t,x(:,2,1))
grid
title('\theta 2 response for \tau 1 = Unit Step')
ylabel('\theta 2(deg)')
xlabel('Time (seconds)')
%return

%
% Tetha1_dot due to \tau 1
%
figure(8)
plot(t,x(:,3,1))
grid
title('theta1dot response for \tau 1 = Unit Step')
ylabel('Tetha1dot (deg/sec)')
% ylabel('$\thta 1 \dot{\gamma} \eta'' (\omega)$', 'Interpreter','latex')
xlabel('Time (seconds)')
%return

%
% Tetha2_dot due to \tau 1
%
figure(9)
plot(t,x(:,4,1))
grid
title('Tetha2dot response for \tau 1 = Unit Step')
ylabel('Tetha2dot (deg/sec)')
xlabel('Time (seconds)')

%return

%
% Step Response due to \tau 2
%

%
% Tetha1 due to \tau 2
%
figure(10)
plot(t,x(:,1,2))
grid
title('\theta 1 response for \tau 2 = Unit Step')
ylabel('\theta 1 (deg)')
xlabel('Time (seconds)')
%return

%
% Tetha2 due to \tau 2
%
figure(11)
plot(t,x(:,2,2))
grid
title('\theta 2 response for Torque2 = Unit Step')
ylabel('\theta 2 (deg))')
xlabel('Time (seconds)')
%return

%
% Theta1_dot due to \tau 2
%
figure(12)
plot(t,x(:,3,2))
grid
title('theta1dot response for \tau 2 = Unit Step')
ylabel('Theta1dot (deg/sec)')
xlabel('Time (seconds)')
%return

%
% Tetha2_dot due to torque 2
%
figure(13)
plot(t,x(:,4,2))
grid
title('theta2dot response for \tau 2 = Unit Step')
ylabel('theta2dot (deg/sec)')
xlabel('Time (seconds)')
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%       Frequency response : Singular values 
%
% u = [T1 (\tau 1) (N-m) T2 (\tau 2) (N-m)]
% x = [theta_1     (deg)     theta_2      (deg)     dot_theta_1 (deg/sec)    dot_theta_2 (deg/sec) ]
% y = [ theta_1     (deg)     theta_2      (deg) ]
%
winit  = -1;
wfin   =  2;
nwpts  = 500;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced fr
%eq points
sv     = sigma(ss(A,B,C,D),w);
sv     = 20*log10(sv);
figure(14); semilogx(w, sv)
%clear sv
title('Plant frequency response Outputs: \theta_1, \theta_2 (deg);      Inputs: \tau_1, \tau_2 (N-m)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
i=14;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Plant SVD Analysis at DC
%
DC = C*inv(-A)*B;
[Udc, Sdc, Vdc ] = svd(DC);
sdc = sigma(ss(DC),w);
sdc = 20*log10(sdc);
i=i+1;
figure(i); semilogx(w, sdc)
title('Plant frequency response at DC Outputs: \theta_1, \theta_2 (deg);      Inputs: \tau_1, \tau_2 (N-m)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)at DC')


% DC =
% 
%    -0.0730    0.0730
%     0.0730   -0.3066
% 
% 
% Udc =
% 
%    -0.2757    0.9612
%     0.9612    0.2757
% 
% 
% Sdc =
% 
%     0.3276         0
%          0    0.0521
% 
% 
% Vdc =
% 
%     0.2757   -0.9612
%    -0.9612   -0.2757

%Maximum singular value of 0.3276 is associated with right singular (input)
%vector Vdc(:,1) and left singular vector Udc(:,1). Since the second
%component of Vdc is much larger than its first component, it follows that
%maximum singular value is primarily associated with second control T2. 
%Since second component of Udc(:,1) is much larger than its first
%component, it follows that maximum singular value is primarily associated
%with the second output \theta 2. The maximum singular value is therefore
%associated with T2 and \theta 2 - the upper link.

%Minimum singular value of 0.0521 is associated with right singular (input)
%vector Vdc(:,2) and left singular vector Udc(:,2). Since the first
%component of Vdc is much larger than its second component, it follows that
%maximum singular value is primarily associated with first control T1. 
%Since first component of Udc(:,1) is much larger than its second
%component, it follows that maximum singular value is primarily associated
%with the first output \theta 1. The minimum singular value is therefore
%associated with T1 and \theta 1 - the lower link.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%
%       Plant SVD analysis at w=5 rad/sec
%
s = 1i*5;
plant_s = C*inv(s*eye(4)-A)*B+D;
[U, S, V] = svd(plant_s);
sw5 = sigma(ss(plant_s),w);
sw5 = 20*log10(sw5);
i=i+1;

figure(i); semilogx(w, sw5)
title('Plant frequency response at \omega = 5rad/s Outputs: \theta_1, \theta_2 (deg);Inputs: \tau_1, \tau_2 (N-m)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)at \omega = 5rad/s')

% 
% plant_s =
% 
%    -0.0086    0.0038
%     0.0038   -0.0190
% 
% 
% U =
% 
%    -0.3113    0.9503
%     0.9503    0.3113
% 
% 
% S =
% 
%     0.0202         0
%          0    0.0074
% 
% 
% V =
% 
%     0.3113   -0.9503
%    -0.9503   -0.3113

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Plant SVD analysis at w=15 rad/sec
s1 = 1i*15;
plant_s1 = C*inv(s1*eye(4)-A)*B+D;
[U1, S1, V1] = svd(plant_s1);
sw15 = sigma(ss(plant_s1),w);
sw15 = 20*log10(sw15);
i=i+1;
figure(i); semilogx(w, sw15)
title('Plant frequency response at \omega = 15rad/s Outputs: \theta_1, \theta_2 (deg);Inputs: \tau_1, \tau_2 (N-m)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)at \omega = 15rad/s')


% plant_s1 =
% 
%    -0.0011    0.0004
%     0.0004   -0.0022
% 
% 
% U1 =
% 
%    -0.3184    0.9480
%     0.9480    0.3184
% 
% 
% S1 =
% 
%     0.0024         0
%          0    0.0009
% 
% 
% V1 =
% 
%     0.3184   -0.9480
%    -0.9480   -0.3184


%SVD analyis at DC, 5rad/sec and 15 rad/sec with increase in frequency we can see 
%more and more coupling effect where T2 effect on \theta 1 increases and T1 effect on \theta 2 increases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%       Control design using LQR method
% Augmenting integrators with plant
ALQ = [ 0*ones(2,2) C
      0*ones(4,2) A ];
  
%   ALQ =
% 
%          0         0    1.0000         0         0         0
%          0         0         0    1.0000         0         0
%          0         0         0         0    1.0000         0
%          0         0         0         0         0    1.0000
%          0         0    4.0421    0.6419         0         0
%          0         0    0.4017    1.7479         0         0

BLQ = [0*ones(2,2)
     B ];
% BLQ =
% 
%          0         0
%          0         0
%          0         0
%          0         0
%     0.2482   -0.0983
%    -0.0983    0.5066
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      LQR design parameter
%Q1=eye(6,6);  %Penalizing all the states
Q1=diag([1 1 1 1 0 0])
%Q1=diag([1 1 1.5 2 0 0])
%Q1=diag([1 1 500 25 0 0]) %Heavy Penalizing Thetha 1 and Thetha 2
R1=0.0001*eye(2,2)
% R1=0.000001*eye(2,2) %Expensive control but good properties
%R1 = 0.01*eye(2,2) %Cheap control but bad properties

[G1, K1, clpoles1] = lqr(ALQ, BLQ, Q1, R1)

% % G1 =
% 
%    99.9999    0.1717  150.4323    8.7327   36.1025    4.9309
%    -0.1717   99.9999    8.2977  125.1246    4.8458   22.9240
% 
% 
% % K1 =
% 
%     1.3245    0.0424    0.3609    0.0489    0.0436    0.0084
%     0.0424    1.2085    0.0491    0.2293    0.0085    0.0214
%     0.3609    0.0491    0.4365    0.0661    0.0663    0.0145
%     0.0489    0.2293    0.0661    0.2578    0.0144    0.0275
%     0.0436    0.0085    0.0663    0.0144    0.0162    0.0041
%     0.0084    0.0214    0.0145    0.0275    0.0041    0.0053
% 
% 
% clpoles1 =
% 
%   -5.2509 + 5.1387i
%   -5.2509 - 5.1387i
%   -3.5620 + 3.0089i
%   -3.5620 - 3.0089i
%   -0.9894 + 0.0000i
%   -0.9999 + 0.0000i

damp(ALQ-BLQ*G1)
                                                                       
%         Pole              Damping       Frequency       Time Constant  
%                                        (rad/TimeUnit)     (TimeUnit)    
%                                                                         
%  -5.25e+00 + 5.14e+00i     7.15e-01       7.35e+00          1.90e-01    
%  -5.25e+00 - 5.14e+00i     7.15e-01       7.35e+00          1.90e-01    
%  -3.56e+00 + 3.01e+00i     7.64e-01       4.66e+00          2.81e-01    
%  -3.56e+00 - 3.01e+00i     7.64e-01       4.66e+00          2.81e-01    
%  -9.89e-01                 1.00e+00       9.89e-01          1.01e+00    
%  -1.00e+00                 1.00e+00       1.00e+00          1.00e+00    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Loop Dynamics 

%Open Loop System
% State x = [z' y' xr']'
% where z is integrator state
%       y is output (tetha 1 and tetha 2)
%       xr contains (\theta 1_dot and \theta 2_dot)

Gz = G1(:,1:2);
Gy = G1(:,3:4);
Gx = G1(:,5:6);

AOL = [0*ones(2,2) 0*ones(2,4)
        -B*Gz  A-B*[0*ones(2,2) Gx]]
    
    
BOL = [-eye(2,2)
        B*Gy ];

COL = [0*ones(2,2) C];

DOL = 0*ones(2,2);
%
%Closed Loop Dynamics
%
ACL = AOL - BOL*COL;
BCL = BOL;
CCL = COL;
DCL = DOL;
cloop = ss(ACL,BCL,CCL,DCL);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Command Following
%
%
t= [0:0.02:5];
[y,t,x] = step(cloop,t);
%
%
%\theta 1 Command
% \theta 1 :r =[1 0] 
%
i=i+1;
figure(i);
plot(t,y(:,:,1))
grid
title('\theta 1 and \theta 2 response to r=[1 0] command')
ylabel('\theta 1 and \theta 2 (deg)')
xlabel('Time (seconds)')

u10 = [-G1 Gy]*[x(:,:,1)'
                ones(1,size(x(:,:,1)')*[0 1]')
                0*ones(1,size(x(:,:,1)')*[0 1]') ];
%Influence on controls due to r=[1 0] Tetha1 command
%
i=i+1;
figure(i);
plot(t,u10)
grid
title('\tau 1 and \tau 2 response to r =[1 0] command')
ylabel('T1, T2 (N-m)')
xlabel('Time (seconds)')
%
% \theta 2 command : r =[0 1] 
%
i=i+1;
figure(i);
plot(t,y(:,:,2))
grid
title('\theta 1 and \theta 2 response to r =[0 1] command')
ylabel('\theta 1 and \theta 2 (deg)')
xlabel('Time (seconds)')
%
%
u01 = [-G1 Gy]*[x(:,:,2)'
                ones(1,size(x(:,:,2)')*[0 1]')
                0*ones(1,size(x(:,:,2)')*[0 1]') ];
%Influence on controls due to r=[0 1] Tetha2 command
%
i=i+1;
figure(i);
plot(t,u01)
grid
title('\tau 1 and \tau 2 response to r =[0 1] command')
ylabel('T1, T2 (N-m)')
xlabel('Time (seconds)')
%
% tetha 1 and \theta 2 comand 
% r =[1 1]
%
[y,t,x] =lsim(cloop,[ones(size(t)) ones(size(t))],t);

i=i+1;
figure(i);
plot(t,y)
grid
title('\theta 1 and \theta 2 response to r=[1 1] Command')
ylabel('\theta 1 and \theta 2 (deg)')
xlabel('TIme (seconds)')
%
% From control : u =[-g gy] [x' xr']'
%
u11= [-G1 Gy]*[x'
               ones(1,size(x')*[0 1]')
               ones(1,size(x')*[0 1]')];
%
i=i+1;
figure(i);
plot(t,u11)
grid
title('\tau 1 and \tau 2 response to r=[1 1]')
ylabel('\tau 1, \tau 2 (N-m)')
xlabel('Time(seconds)')
%
% r=[-1 1] command
%
[y,t,x] =lsim(cloop,[-ones(size(t)) ones(size(t))],t);
i=i+1;
figure(i)
plot(t,y)
grid
title('\theta 1 and \theta 2 response to r=[-1 1] command')
ylabel('\theta 1 and \theta 2 (deg)')
xlabel('Time (seconds)')
%
% Form control :u=[-g gy][x' xr']'
%
um11 = [-G1 Gy]*[x'
                -ones(1,size(x')*[0 1]')
                 ones(1,size(x')*[0 1]')];
i=i+1;
figure(i)
plot(t,um11)
grid
title('\tau 1 and \tau 2 response to r=[-1 1] Command')
ylabel('\tau 1, \tau 2 (N-m)')
xlabel('Time (seconds)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% LQ OPEN LOOP FREQUENCY RESPONSE 
%
w = logspace(-1,2,500);
sv = sigma(ss(ALQ, BLQ, G1, 0*ones(2,2)),w);
sv = 20*log10(sv);
i=i+1;
figure(i);semilogx(w, sv)
%clear sv
title('Open Loop Singular Values: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LQ CLOSED LOOP FREQUENCY RESPONSE 
%
w = logspace(-1,2,500);
sv = sigma(ss(ALQ-BLQ*G1, BLQ, -G1, eye(2,2)-0*ones(2,2)),w);
sv = 20*log10(sv);
i=i+1;
figure(i);
semilogx(w, sv)
%clear sv
title('LQ  Sensitivity: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
%return

i=i+1;
figure(i)
sv =sigma(ss(ACL,BCL,-CCL,eye(2,2)-DCL),w);
sv=20*log10(sv);
semilogx(w,sv)
title('LQ Sensitivity: Error Signal')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

sv = sigma(ss(ACL,BCL,CCL,DCL),w);
sv = 20*log10(sv);
i=i+1;
figure(i);
semilogx(w, sv)
title('LQ Complementary Sensitivity: Plant Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


sv = sigma(ss(ALQ-BLQ*G1, BLQ, G1, 0*ones(2,2)),w);
sv = 20*log10(sv);
i=i+1;
figure(i);
semilogx(w, sv)
%clear sv
title('LQ Complementary Sensitivity: Plant Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%% Every procedure same till line 602
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Dr. Rodriguez 's LQR code
%       Control using LQG method, Loop shaping
%Augmenting plant with integrators to achieve zero steady state error

a= [zeros(nc,nc+ns); B A];
b= [eye(nc,nc);zeros(ns,nc)];
c= [zeros(nc,nc) C];
d= [zeros(2,2)];
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(a, b, c, d),w);
sv     = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Design Plant Singular Values (LQG)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
%
%% Using Kalman Filtering Idea to design target loop (at output)
%   Design1
Ll = -inv(C*inv(A)*B);
Lh = -inv(A)*B*Ll;
L=[Ll;Lh];
z=1.2;



winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
gfolsv = sigma(ss(a, L, c, d),w);
gfolsv = 20*log10(gfolsv);
i=i+1;
figure(i); semilogx(w, gfolsv)
%clear sv
title('G_{FOL} Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mu = 0.1; %Target Loop bandwidth parameter
%Decreasing mu will raise the bandwidth
%Increasing mu will decrease the bandwidth
%mu=0.01;
%mu=10;

x   = are(a',(1/mu)*c'*c,L*L')  %Algebraic riccati equation
h   = (1/mu)*x*c';                % Kalman gain matrix H_f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%   Target closed loop poles and zeros
tpoles = eig(a-h*c)
tzeros = tzero(a-h*c,h,c,d)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%   Pre-Filter
%
transf = tf({z 0;0 z}, {[1 z] 1; 1 [1 z]});
fil =ss(transf);
%
%   Target Open Loop SVD at s=j20
%
s =j*20
TOL_s = c*inv(s*eye(size(a))-a)*h+d;
[TOLu, TOLs, TOLv] = svd(TOL_s)

%
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(a,h,c,d),w);
tlsv   = 20*log10(sv);
i=i+1;
figure(i);
semilogx(w, tlsv)
title('Target Open Loop Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Target Sensitivity Singular Values
%
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(a-h*c,h,-c,eye(2,2)),w);
sv     = 20*log10(sv);
i=i+1; figure(i)
semilogx(w, sv)
title('Target Sensitivity Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Target Complementary Sensitivity Singular Values
%
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(a-h*c,h,c,d),w);
sv     = 20*log10(sv);
i=i+1; figure(i)
semilogx(w, sv)
title('Target Comp Sensitivity Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Target Closed Loop Singular Values
%
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
[ntacl, ntbcl, ntccl, ntdcl ] = series(fil.a, fil.b, fil.c, fil.d, a-h*c,h,c,d);
sv     = sigma(ss(ntacl, ntbcl, ntccl, ntdcl),w);
sv     = 20*log10(sv);
i=i+1; figure(i)
semilogx(w, sv)
title('Target Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('T_{ry}, Singular Values (dB)')

%
% Target Step Responses
%
tinit = 0;
tinc  = 0.005;
tfin  = 4.0;
t     = [tinit:tinc:tfin]';                   % Vector of uniformly spaced time points
r1    = [ones(size(t)) zeros(size(t))];
r2    = [zeros(size(t)) ones(size(t))];
ty1   = lsim(ss(ntacl, ntbcl, ntccl, ntdcl), r1,t);
ty2   = lsim(ss(ntacl, ntbcl, ntccl, ntdcl), r2,t);
i=i+1; figure(i)
plot(t,ty1(:,1),'b', t,ty1(:,2),'r')
grid
title('Target Responses to \theta_1 Step Reference Command')
ylabel('\theta_1, \theta_2 (degrees)')
xlabel('Time (seconds)')

i=i+1; figure(i)
plot(t,ty2(:,1),'b', t,ty2(:,2),'r')
grid
title('Target Responses to \theta_2 Step Reference Command')
ylabel('\theta_1, \theta_2 (degrees)')
xlabel('Time (seconds)')

%%
% Loop Transfer Recovery at Output
%
q            = c'*c;
rho          = 1e-13;	                      % LQG/LTRO Recovery Parameter
[g, kp, clp] = lqr(a, b, q, rho*eye(nc,nc));  % Compute Control Gain Matrix G


ak = [0*ones(nc,2)  g; zeros(ns+nc,2)  a-b*g-h*c];
bk = [zeros(nc,2); h];
ck = [eye(nc,2)  zeros(nc, ns+nc)];
dk = [zeros(2,2)];
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
kpoles = eig(ak)                               % Compensator Poles
kzeros = tzero(ak, bk, ck, dk)                 % Compensator Zeros

winit  = -2;
wfin   = 4;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);
sv     = sigma(ss(ak, bk, ck, dk),w);
sv     = 20*log10(sv);
i=i+1; figure(i); semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Form Open Loop System
%
[al, bl, cl, dl ] = series(ak, bk, ck, dk, A, B, C, D);
  
olpoles = eig(al)                          % Open Loop Poles
olzeros = tzero(al,bl,cl,dl)               % Open Loop Zeros
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(al, bl, cl, dl), w);
sv      = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv, w, tlsv)
%clear sv
title('Open Loop Singular Values at Error (Recovered and Target)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Form Open Loop System
%
[ali, bli, cli, dli ] = series(A,B,C,D, ak, bk, ck, dk);
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(ali, bli, cli, dli ), w);
sv      = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Open Loop Singular Values at Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Form Closed Loop System
%
acl     = al-bl*cl;
bcl     = bl;
ccl     = cl;
dcl     = dl;

clpoles = eig(acl)               % Closed Loop Poles
damp(clpoles)
clzeros = tzero(acl,bcl,ccl,dcl) % Closed Loop Zeros (r to y)

%
% Sensitivity at Error
%
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(acl, bcl, -ccl, eye(2)),w);
sv      = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values at Error')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Sensitivity at Input
%
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(ali-bli*cli, bli, -cli, eye(2)),w);
sv      = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values at Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Complementary Sensitivity
%
winit   = -1;
wfin    = 2;
nwpts   = 200;
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(acl, bcl, ccl, dcl),w);
sv      = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values at Ouput')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Add Pre-filter to Closed Loop Transfer Function Matrix From r to y
%

[nacl, nbcl, nccl, ndcl ]   = series(fil.a, fil.b, fil.c, fil.d, acl,bcl,ccl,dcl);
sv                          = sigma(ss(nacl, nbcl, nccl, ndcl),w);
sv                          = 20*log10(sv);

i=i+1;figure(i); semilogx(w, sv)
%clear sv
title('Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Form Reference to Controls (r to u)
%
[aru, bru, cru, dru ] = series(acl, bcl, -ccl, eye(2), ak, bk, ck, dk);
[aru, bru, cru, dru ] = series(fil.a, fil.b, fil.c, fil.d, aru, bru, cru, dru );    
winit                 = -1;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(aru, bru, cru, dru),w);
sv                    = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Reference to Control Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Form Input Disturbance to Output
%
[ady, bdy, cdy, ddy ] = series(A,B,C,D, acl, bcl, -ccl, eye(2));
winit                 = -2;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(ady, bdy, cdy, ddy),w);
sv                    = 20*log10(sv);
i=i+1;
figure(i); semilogx(w, sv)
%clear sv
title('Input Disturbance to Output Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%
% Input Disturbance Analysis at DC
%
s = j*7
tdy              = cdy*inv(s*eye(size(ady))-ady)*bdy +  ddy;
[utdy stdy vtdy] = svd(tdy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%   Closed Loop Step Responses
%
y1    = lsim(ss(nacl,nbcl, nccl, ndcl), r1, t); 
i=i+1;
figure(i);
plot(t,y1, t, ty1)
plot(t,y1)
grid
title('Response to Step Reference Command for \theta_1')
ylabel('Outputs (deg)')
xlabel('Time (seconds)')

y2    = lsim(ss(nacl,nbcl, nccl, ndcl), r2, t); 
i=i+1;
figure(i); plot(t,y2,  t, ty2)
plot(t,y2)
grid
title('Response to Step Reference Command for \theta_2')
ylabel('Outputs (deg)')
xlabel('Time (seconds)')

u1    = lsim(ss(aru, bru, cru, dru), r1, t); 
i=i+1;
figure(i); plot(t,u1)
grid
title('Response to Step Reference Command for \theta_1')
ylabel('Controls (N-m)')
xlabel('Time (seconds)')

u2    = lsim(ss(aru, bru, cru, dru), r2, t); 
i=i+1;
figure(i);plot(t,u2)
grid
title('Response to Step Reference Command for \theta_2')
ylabel('Controls (N-m)')
xlabel('Time (seconds)')





