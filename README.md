clear;clc

X=xlsread('D:\文件\强降雨作用下破坏概率演化\案例\龙家台滑坡\龙家台.xlsx','A2:L41');

m1=10000;
n=40;

q=4/3.6*10^(-5);
theta_delta=normrnd(0.19,0.0019,m1,n);
gama_s=20.5;
gama_a=18.6;
phi_b=0.226892802759263;
phi=lognrnd(-1.34022676143877,0.00999975001354021,m1,n);
c=lognrnd(2.990757108127,0.00999975001354021,m1,n);
s=lognrnd(4.78744174528188,0.00999975001354021,m1,n);
psai=lognrnd(-0.223193548814376,0.00999975001354021,m1,n);
Ks=lognrnd(-12.7169482667963,0.00999975001354021,m1,n);
gama_w=10;
t=0;


for j=1:n
H(j)=X(j,2)-X(j,8);
h(j)=X(j,10)-X(j,8);
l(j)=((X(j,7)-X(j,5))^2+(X(j,8)-X(j,6))^2)^0.5;
alpha(j)=atan((X(j,6)-X(j,8))/(X(j,5)-X(j,7)));
beta(j)=atan((X(j,4)-X(j,2))/(X(j,3)-X(j,1)));
end

alpha(n+1)=0;
h(n+1)=0;
H(n+1)=0;
El(1)=0;
Eu(n+1)=0;


for i=1:m1

for j=1:n

z_p(i,j)=psai(i,j)/(q/Ks(i,j)-1)/cos(beta(j));
t_p(i,j)=psai(i,j)*theta_delta(i,j)/((q/Ks(i,j)-1)*q*(cos(beta(j)))^2);

if t<t_p
z_f(i,j)=q*cos(beta(j))*t/theta_delta(i,j);
else
syms z
eqn=theta_delta(i,j)/Ks(i,j)/cos(beta(j))*(z-z_p(i,j)-psai(i,j)/cos(beta(j))*log((z*cos(beta(j))+psai(i,j))/(z_p(i,j)*cos(beta(j))+psai(i,j))))+t_p(i,j)-t;
zf=double(solve(eqn,z));
z_f(i,j)=max(zf);
end

if 1<j<n
eta(i,j)=min((z_f(i,j)/cos(beta(j))+0.5*(h(j)+h(j+1)))/(0.5*(H(j)+H(j+1))),1);
elseif j==1
eta(i,j)=min(1-(1-z_f(i,j)/cos(beta(j))/H(j+1))^2+h(j+1)/H(j+1),1);
else
eta(i,j)=min(1-(1-z_f(i,j)/cos(beta(j))/H(j))^2+h(j)/H(j),1);
end


c_psai(i,j)=c(i,j)+s(i,j)*tan(phi_b)*(1-eta(i,j))^3.4;

G(j)=0.5*(gama_s-gama_w)*(2*z_f(i,j)/cos(beta(j))+h(j)+h(j+1))*l(j)*cos(alpha(j))+0.5*gama_a*(H(j)+H(j+1)-2*z_f(i,j)/cos(beta(j))-h(j)-h(j+1))*l(j)*cos(alpha(j));
P(j)=0.5*gama_w*l(j)*(h(j)+h(j+1))*sin(alpha(j))*cos(alpha(j));
Q(j)=gama_w*l(j)*z_f(i,j)*tan(beta(j))*sin(alpha(j));

lama(j)=cos(alpha(j)-alpha(j+1))-sin(alpha(j)-alpha(j+1))*tan(phi(i,j));
epsi(j)=cos(alpha(j))*tan(phi(i,j))-sin(alpha(j));
kapa(j)=cos(alpha(j)-beta(j))+sin(alpha(j)-beta(j))*tan(phi(i,j));

end


for k=2:n
El(k)=1/lama(k-1)*(El(k-1)+epsi(k-1)*G(k-1)-kapa(k-1)*Q(k-1)-P(k-1)+c_psai(i,k-1)*l(k-1));
if El(k)<0
El(k)=0;
end
end

for k=(n-1):(-1):1
Eu(k+1)=lama(k+1)*Eu(k+2)-epsi(k+1)*G(k+1)+kapa(k+1)*Q(k+1)+P(k+1)-c_psai(k+1)*l(k+1);
if Eu(k+1)<0
Eu(k+1)=0;
end
end

for j=1:1:n
D(i,j)=(Eu(j+1)*cos(alpha(j)-alpha(j+1))-El(j)+G(j)*sin(alpha(j))+P(j)+Q(j)*cos(alpha(j)-beta(j)));
R(i,j)=((G(j)*cos(alpha(j))-Eu(j+1)*sin(alpha(j)-alpha(j+1))-Q(j)*sin(alpha(j)-beta(j)))*tan(phi(i,j))+c_psai(i,j)*l(j));
FLSM(i,j)=(R(i,j)-D(i,j))/R(j);
end

end
