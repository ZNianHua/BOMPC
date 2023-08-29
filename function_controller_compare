function [a,omega,kappa_i_1] = controller_compare(state_i,sequence_i_1,d0,Tstep,kappa_last,K)
%恒间距策略，论文中h_i趋于0
if nargin<6
    k1=3;
    k2=3;
else
    k1=K(1);
    k2=K(2);
end
h_i=1e-3;
x_i=state_i(1);
y_i=state_i(2);
theta_i=state_i(3);
v_i=state_i(4);
x_i_1=sequence_i_1(1,end);
y_i_1=sequence_i_1(2,end);
theta_i_1=sequence_i_1(3,end);
v_i_1=sequence_i_1(4,end);
omega_i_1=sequence_i_1(6,end);

if v_i_1==0
    kappa_i_1=0;
else
    kappa_i_1=abs(omega_i_1/v_i_1);
end
kappa_i_1_dot=(kappa_i_1-kappa_last)/Tstep;
if kappa_i_1==0%论文参数overwilde s_i-1——式18
    over_s_i_1=0;
else
    over_s_i_1=(-1+sqrt(1+kappa_i_1^2*(d0+h_i*v_i)))/kappa_i_1;
end

alpha_i=atan(kappa_i_1*(d0+h_i*v_i));%论文α_i——式20
s_ki=(1-cos(alpha_i))/(kappa_i_1^2);%论文参数s_ki——式26
s_ai=h_i*sin(alpha_i);%论文s_ai——式26
%论文矩阵β_i1——式33
beta_1i=[-sin(theta_i);cos(theta_i)]*v_i*tan(alpha_i)+[cos(theta_i),-sin(theta_i);sin(theta_i),cos(theta_i)]*[over_s_i_1*omega_i_1;-s_ki*kappa_i_1_dot];
%论文参数μ_i——式37
miu_i=h_i*(d0+h_i*v_i)*(1-sin(alpha_i)*sin(theta_i_1-theta_i));
%论文矩阵Kappa_12i的逆——式36
Kappa_12i=[(d0+h_i*v_i)*cos(theta_i),(d0+h_i*v_i)*sin(theta_i);-h_i*sin(theta_i)-s_ai*cos(theta_i_1),h_i*cos(theta_i)-s_ai*sin(theta_i_1)]./miu_i;
%论文z_1i~z_4i——式19
z_1i=x_i_1+over_s_i_1*sin(theta_i_1)-x_i-(d0+h_i*v_i)*cos(theta_i);
z_2i=y_i_1-over_s_i_1*cos(theta_i_1)-y_i-(d0+h_i*v_i)*sin(theta_i);
z_3i=v_i_1*cos(theta_i_1)-v_i*cos(theta_i+alpha_i);
z_4i=v_i_1*sin(theta_i_1)-v_i*sin(theta_i+alpha_i);
%控制量计算——式35
control=Kappa_12i*([k1*z_1i;k2*z_2i]+[z_3i;z_4i]./(cos(alpha_i))+beta_1i);
a=control(1);
omega=control(2);
%%控制量限制

end
