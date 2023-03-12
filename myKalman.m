function [X_hatrec, Y_hatrec, K_rec, P_rec] = myKalman(A,B,C,D,x0,u,y,P0,Qo,Ro)
Size = length(u);
X_hatrec = zeros(length(A), Size);
Y_hatrec = zeros(length(D), Size);
P_rec = zeros(Size, 2*length(A));
K_rec = zeros(Size, length(A));
X_hat = x0;
Pp = P0;
for i = 1:Size
 X_hatrec(:,i) = X_hat;
 P_rec(i,:) = [Pp(1,:), Pp(2,:)];
 K = A * Pp * C' / (C * Pp * C' + Ro);
 K_rec(i,:) = [K(1,:), K(2,:)];
 X_hat = (A - K * C) * X_hat + B * u(i) + K * y(:,i);
 Pp = A * Pp * A' - A * Pp * C' / (C * Pp * C' + Ro) * C * Pp * A'
 + Qo;
 Y_hatrec(:,i) = (C * X_hat);
end
end
