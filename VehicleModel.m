function [state_next,Power_vector] = VehicleModel(state_pre,vehicle_model_param,U,phi,Tstep,uncertain_flag)
if nargin<6
    uncertain_flag=zeros(7,1);
end

M=vehicle_model_param(1);
Iz=vehicle_model_param(2);
a=vehicle_model_param(3);
b=vehicle_model_param(4);
Cx=vehicle_model_param(5);
Cy=vehicle_model_param(6);
Cf=vehicle_model_param(7);
Cr=vehicle_model_param(8);
hg=vehicle_model_param(9);
f=vehicle_model_param(10);
tao=vehicle_model_param(11);
miu_f=vehicle_model_param(12);
miu_r=vehicle_model_param(13);
R=vehicle_model_param(14);
i_0=vehicle_model_param(15);
yita=vehicle_model_param(16);
g=vehicle_model_param(17);
 

L=a+b;

delta=U(2);
T_t=U(1);

Xp=state_pre(1);
Yp=state_pre(2);
thetap=state_pre(3);
Vxp=state_pre(4);
Vyp=state_pre(5);
theta_dotp=state_pre(6);
accp=state_pre(7);

Fzf=M*g*b*cos(phi)/L-M*g*hg*sin(phi)/L-M*hg*accp/L;
Fzr=M*g*a*cos(phi)/L+M*g*hg*sin(phi)/L+M*hg*accp/L;


if Vxp==0||theta_dotp==0
    alpha_f=0;
    alpha_r=0;
else
    alpha_f=-(delta-atan2((Vyp+a*theta_dotp),Vxp));
    alpha_r=(atan2((Vyp-b*theta_dotp),Vxp));
end
Fyf=Cf*alpha_f;
Fyr=Cr*alpha_r;
Fxf=T_t*yita*i_0/R*miu_f/(miu_r+miu_f);
Fxr=T_t*yita*i_0/R*miu_r/(miu_r+miu_f);
F_air=Cx*abs(Vxp)*Vxp;
F_roll=M*g*f*cos(delta);
Fx=Fxr+Fxf*cos(delta)-Fyf*sin(delta)-M*g*sin(phi)-F_roll-F_air;
Fy=Fxf*sin(delta)+Fyf*cos(delta)+Fyr;




Vx_dot=Fx/M;
Vy_dot=Fy/M;
theta_dot_dot=(a*Fxf*sin(delta)+a*Fyf*cos(delta)-b*Fyr)/Iz;
acc_dot=(Vx_dot-accp)/(tao);

X=Xp+Tstep*Vxp*cos(thetap)-Tstep*Vyp*sin(thetap)+uncertain_flag(1)*(rand(1)-0.5)*1;
Y=Yp+Tstep*Vxp*sin(thetap)+Tstep*Vyp*cos(thetap)+uncertain_flag(2)*(rand(1)-0.5)*1;
theta=thetap+Tstep*theta_dotp+uncertain_flag(3)*rand(1)*0.1;
Vx=Vxp+Tstep*Vx_dot-Vyp*theta_dotp*Tstep+uncertain_flag(4)*(rand(1)-0.5)*0.5;
Vy=Vyp+Tstep*Vy_dot+Vxp*theta_dotp*Tstep+uncertain_flag(5)*(rand(1)-0.5)*0.1;
theta_dot=theta_dotp+Tstep*theta_dot_dot+uncertain_flag(6)*(rand(1)-0.5)*0.1;
acc=accp+acc_dot*Tstep+uncertain_flag(7)*(rand(1)-0.5)*0.1;
% acc=Vx_dot+uncertain_flag(7)*(rand(1)-0.5)*0.1;
P_x=M*accp*Vxp;
P_y=M*Vy_dot*Vyp;
P_roll=F_roll*Vxp;
P_air=F_air*Vxp;
if theta>pi
    theta=theta-2*pi;
elseif theta<-pi
    theta=theta+2*pi;
end

state_next=[X;Y;theta;Vx;Vy;theta_dot;acc];
Power_vector=[P_x,P_y,P_roll,P_air]';
end

