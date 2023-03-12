function [X,Y]=myss(A,B,C,D,G,H,u,x0,v,d)
Size=length(u);
X=zeros(length(A),Size);
Y=zeros(length(D),Size);
X(1:length(A),1) = x0;
for i = 1:Size - 1
 Y(1:length(D),i) = C * X(1:length(A), i) + D * u(i) + H * v(:,i);
 X(1:length(A), i+1) = A * X(1:length(A), i) + B * u(i) + G *
 d(:,i);
end
end