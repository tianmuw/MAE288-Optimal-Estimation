clear;clc;
%% Wind Model
[numb, denb]=butter(2,.05);
t=[0:.5:200]';
[T,~] = size(t);
Qs_det=1.3;
u_w = randn(T,1)*Qs_det;
wi=filter(numb,denb,u_w); % wi is the output of wind velocity.
numb=numb(1:2); % This is the delay process, making the Ds zero.
[As,Bs,Cs,Ds]=tf2ss(numb,denb); % The wind system As,Bs,Cs,Ds.
wind_force_ratio = 10; %This is a coefficient of relation between wind velocity and force.
%% Cart Model
[Ac,Bc,Cc,Dc]=ssdata(c2d(ss([0 1;0 -.2],[0;1],eye(2),zeros(2,1)),.5)); % This is the cart model given by the prof.
Gc=[0;1];
Motor = sin(2*pi*t/10);
%% Combined Model
F = [As/wind_force_ratio, zeros(2); Gc*Cs/Wind_force_ratio, Ac];
G = [zeros(2,1);Bc];
B1 = [Bs zeros(2,3);zeros(2,3) Gc*Ds]/Wind_force_ratio;
H = [zeros(2), Cc];
J = [0;0];
B2 = eye(2);
%% Compute The Real Data
X0 = [0.5;0.5;2;1];
% The measurement noise v is a 2x1 vector, the first element represents the
% postion and the second element represents the velocity.
% v = [v_pos; v_vel];
% To get a v.
v_pos = 0.05;
v_vel = 0.05;
v_all = [v_pos; v_vel];
% To realize these measurement noise.
v_pos_real = randn(1,T)*sqrt(v_pos);
v_vel_real = randn(1,T)*sqrt(v_vel);
v_all_real = [v_pos_real; v_vel_real];

% The process noise is actually the noise in wind model, which is a randn
% number.
Qs=3;
wind_process_noise = randn(2,T)*Qs;
wind_meas_noise = randn(2,T);
d = [wind_process_noise; wind_meas_noise];
Ro = diag([v_pos,v_vel]);
Qs = 5;
Qo = diag([Qs,0,1e-5,1e-5]);
[L,~,~] = dlqr(F,G,Qo,v_vel);
Fnew = F - G*L;
Gnew = zeros(4,1);
Hnew = H;
Jnew = J;
[X_ss, Y_ss] = myss(Fnew,Gnew,Hnew,Jnew,B1,B2,Motor,X0,v_all_real,d);

%% The Kalman Filter Process with Feedback-Control
% We need to get Ro and Qo.
% Ro is the measurement noise matrix which is 2x2.
% Qo is the process noise matrix which is 4x4;
% In the combined model, the measurement noise only has v_pos and v_vel,
% thus the Ro should be diag([v_pos v_vel])
% Given that the process noise only has the wind noise Qs. Thus Qo should
% be diag([0,Qs,0,0]);
% Get X0 and P0
x0 = [0.5;0.5;2;1];
P0 = eye(4);
[X_hatrec, Y_hatrec, K_rec, P_rec] = myKalman(Fnew,Gnew,Hnew,Jnew,x0,Motor,Y_ss,P0,Qo,Ro);
%% Plot
figure(1);
plot(t,Y_hatrec(1,:));
hold on;
plot(t,Y_ss(1,:));
hold off;
legend('Kalman','Real data');
title('The position data');
figure(2);
plot(t,Y_hatrec(2,:));
hold on;
plot(t,Y_ss(2,:));
hold off;
legend('Kalman','Real data');
title('The velocity data');
Y_tilda = Y_hatrec - Y_ss;
figure(3);
plot(xcorr(Y_tilda(1,:)));
title('The covariance of Y in position');
figure(4);
plot(xcorr(Y_tilda(2,:)));
title('The covariance of Y in velocity');