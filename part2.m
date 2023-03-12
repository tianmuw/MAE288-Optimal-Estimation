clear;clc;
%% Wind Model
[numb, denb]=butter(2,.05);
t=[0:.5:200]';
[T,~] = size(t);
Qs_det=1e-3;
u_w = randn(T,1)*Qs_det;
wi=filter(numb,denb,u_w); % wi is the output of wind velocity.
numb=numb(1:2); % This is the delay process, making the Ds zero.
[As,Bs,Cs,Ds]=tf2ss(numb,denb); % The wind system As,Bs,Cs,Ds.
Wind_force_ratio = 10; %This is a coefficient of relation between wind velocity and force.

%% Cart Model
[Ac,Bc,Cc,Dc]=ssdata(c2d(ss([0 1;0 -.2],[0;1],eye(2),zeros(2,1)),.5)); % This is the cart model given by the Prof.
Gc=[1e-2;10];
Motor = sin(2*pi*t/10);

%% Combined Model
F = [As/Wind_force_ratio, zeros(2); Gc*Cs/Wind_force_ratio, Ac];
G = [zeros(2,1);Bc];
B1 = [Bs zeros(2,3);zeros(2,3) Gc*Ds]/Wind_force_ratio;
H = [zeros(2), Cc];
J = [0;0];
B2 = eye(2);

%% Compute The Real Data
X0 = [0;0;2;1];
% The measurement noise v is a 2x1 vector, the first element represents the
% postion and the second element represents the velocity.
% v = [v_pos; v_vel];
% To get a v.
v_pos = 0.5;
v_vel = 0.5;
v_all = [v_pos; v_vel];
% To realize these measurement noise.
v_pos_real = randn(1,T)*sqrt(v_pos);
v_vel_real = randn(1,T)*sqrt(v_vel);
v_all_real = [v_pos_real; v_vel_real];

% The process noise is actually the noise in wind model, which is a randn
% number.
Qs=1e-3;
wind_process_noise = randn(T,2)*Gc*sqrt(Qs);
wind_meas_noise = randn(2,T);
d = [wind_process_noise';zeros(1,T); wind_meas_noise];
[X_ss, Y_ss] = myss(F,G,H,J,B1,B2,Motor,X0,v_all_real,d);

%% The Kalman Filter Process
% We need to get Ro and Qo.
% Ro is the measurement noise matrix which is 2x2.
% Qo is the process noise matrix which is 4x4;
% In the combined model, the measurement noise only has v_pos and v_vel,
% thus the Ro should be diag([v_pos v_vel])
% Given that the process noise only has the wind noise Qs. Thus Qo should
% be diag([0,Qs,0,0]);
Ro = diag([v_pos, v_vel]);
Qs = 1e-2;
Qo = diag([7,Qs,1e-4,1e-5]);
% Get X0 and P0
x0 = [0;0;2;1];
P0 = eye(4);
[X_hatrec, Y_hatrec, K_rec, P_rec] = myKalman(F,G,H,J,x0,Motor,Y_ss,P0,Qo,Ro);
%% Plot
figure(1);
plot(t, Y_ss(1,:));
hold on
plot(t, Y_hatrec(1,:));
plot(t, Y_ss(2,:))
plot(t, Y_hatrec(2,:));
hold off;
legend('The real data of position','Kalman data of position','The real data of velocity','Kalman data of velocity');
title('The position and velocity data');
% figure(2);
% plot(t, X_ss(3,:), t, X_hatrec(3,:));
% hold on;
% plot(t, X_ss(4,:), t, X_hatrec(4,:));
%Y_tilda = Y_ss - Y_hatrec;
Y_tilda = Y_hatrec - Y_ss;
%X_tilda = X_hatrec - X_ss;
%Y_tilda_2_c=xcorr(Y_tilda(2,:)',Y_tilda(2,:)');
figure(3);
plot(xcorr(Y_tilda(1,:)));
title('Covariance of Y in position');
figure(4);
plot(xcorr(Y_tilda(2,:)));
title('Covariance of Y in velocity');
X_tilda = X_hatrec - X_ss;
figure(5);
plot(xcorr(X_tilda(3,:)));
title('Covariance of X in position');
figure(6);
plot(xcorr(X_tilda(4,:)));
title('Covariance of X in velocity');