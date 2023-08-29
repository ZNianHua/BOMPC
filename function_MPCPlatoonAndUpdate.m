function [state_next,control_applied,MPC_result,Power_vector] = MPCPlatoonAndUpdate(state_now,traj_ref,vehicle_model_param,phi,Tstep,MPC_hori,weight_vector,MPC_U_input_total_last)
%% 数据处理
Np=MPC_hori(1);
Nc=MPC_hori(2);
[Nx,~]=size(traj_ref);
% Kd_X=weight_vector(1);
% Kd_Y=weight_vector(2);
% Ktheta=weight_vector(3);
% Kv_x=weight_vector(4);
% Kv_y=weight_vector(5);
% Komega=weight_vector(6);
% Ka=weight_vector(7);
% K_Tq=weight_vector(8);
% K_delta=weight_vector(9);
Q_weight=diag(weight_vector(1:Nx));
R_weight=diag(weight_vector(Nx+1:Nx+2));
E_weight=weight_vector(end);
%% 车辆模型
import casadi.*
Nu=2;N_power=4;
X_state = MX.sym('x', Nx);
U_input = MX.sym('u', Nu);
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
delta=U_input(2);
T_t=U_input(1);
Xp=X_state(1);
Yp=X_state(2);
thetap=X_state(3);
Vxp=X_state(4);
Vyp=X_state(5);
theta_dotp=X_state(6);
accp=X_state(7);
Fzf=M*g*b*cos(phi)/L-M*g*hg*sin(phi)/L-M*hg*accp/L;
Fzr=M*g*a*cos(phi)/L+M*g*hg*sin(phi)/L+M*hg*accp/L;
alpha_f=-(delta-atan2((Vyp+a*theta_dotp),Vxp));
alpha_r=(atan2((Vyp-b*theta_dotp),Vxp));
Fyf=Cf*alpha_f;
Fyr=Cr*alpha_r;
Fxf=T_t*yita*i_0/R*miu_f/(miu_r+miu_f);
Fxr=T_t*yita*i_0/R*miu_r/(miu_r+miu_f);
F_air=Cx*abs(Vxp)*Vxp;
F_roll=M*g*f*cos(delta);
Fx=Fxr+Fxf*cos(delta)-Fyf*sin(delta)-M*g*sin(phi)-F_roll-F_air;
Fy=Fxf*sin(delta)+Fyf*cos(delta)+Fyr;
Vx_dot=Fx/M+Vyp*theta_dotp;
Vy_dot=(Fy/M-Vxp*theta_dotp);
theta_dot_dot=(a*Fxf*sin(delta)+a*Fyf*cos(delta)-b*Fyr)/Iz;
acc_dot=(Vx_dot-accp)/(tao);
X=Xp+Tstep*Vxp*cos(thetap)-Tstep*Vyp*sin(thetap);
Y=Yp+Tstep*Vxp*sin(thetap)+Tstep*Vyp*cos(thetap);
theta=thetap+Tstep*theta_dotp;
Vx=Vxp+Tstep*accp-Vyp*theta_dotp*Tstep;
Vy=Vyp+Tstep*Vy_dot+Vxp*theta_dotp*Tstep;
theta_dot=theta_dotp+Tstep*theta_dot_dot;
acc=accp+acc_dot*Tstep;
P_x=M*accp*Vxp;
P_y=M*Vy_dot*Vyp;
P_roll=F_roll*Vxp;
P_air=F_air*Vxp;
Power_vector=vertcat(P_x,P_y,P_roll,P_air);
X_next=vertcat(X, Y, theta, Vx, Vy, theta_dot, acc);
f_state= Function('f', {X_state, U_input}, {X_next});
f_power=Function('f_power', {X_state, U_input}, {Power_vector});
l_X_state = MX.sym('x', Nx);
l_X_ref = MX.sym('x', Nx);
l_U_input = MX.sym('u', Nu);
l_U_last = MX.sym('u', Nu);
l_state=(l_X_state-l_X_ref)'*Q_weight*(l_X_state-l_X_ref);
l_delta_U=(l_U_input-l_U_last)'*R_weight*(l_U_input-l_U_last);
l_power_vector=f_power(l_X_state,l_U_input);
l_power=l_power_vector'*eye(N_power)*l_power_vector;
loss_step=l_state+l_delta_U+l_power*E_weight;
% loss_step=l_state+l_delta_U+0;
l=Function('l',{l_X_state,l_X_ref,l_U_input,l_U_last}, {loss_step});
%% MPC 问题
% MPC_U_input_total_last=zeros(Nu*Nc,1);
% MPC_U_input_total_last(1:Nu)=Control_last;
MPC_X_state = cell(Np,1);
MPC_X_state{1}=state_now;
MPC_problem=casadi.Opti();
MPC_U_input_total=MPC_problem.variable(Nu*Nc);
MPC_problem.set_initial(MPC_U_input_total,MPC_U_input_total_last)
J=l(state_now,traj_ref(:,1),0,0);
l_U_last=MPC_U_input_total_last(1:Nu);
for j=1:Np-1
    if j<Nc+1
%         MPC_U_input=MPC_problem.parameter('MPC_U_input',Nu);
        MPC_U_input=MPC_U_input_total(Nu*j-Nu+1:Nu*j);
        MPC_X_state{j+1}=f_state(MPC_X_state{j},MPC_U_input);
        J=J+l(MPC_X_state{j+1},traj_ref(:,j+1),MPC_U_input,l_U_last);
        l_U_last=MPC_U_input;
    else
        MPC_U_input=l_U_last;
        MPC_X_state{j+1}=f_state(MPC_X_state{j},MPC_U_input);
        J=J+l(MPC_X_state{j+1},traj_ref(:,j+1),MPC_U_input,l_U_last);
    end
end

T_power_max=150;
steer_max=35;
delta_steer_max=(steer_max/180)*pi;
Tq_max=T_power_max*i_0*yita;
lb_1=zeros(1,Nc)-Tq_max;
lb_2=zeros(1,Nc)-delta_steer_max;
ub_1=zeros(1,Nc)+Tq_max;
ub_2=zeros(1,Nc)+delta_steer_max;
lb=[lb_1;lb_2];
lb=lb(:);
ub=[ub_1;ub_2];
ub=ub(:);
MPC_problem.subject_to(lb<=MPC_U_input_total<=ub);
MPC_problem.minimize(J);
plugin_opts = struct('expand',true,'print_time',0);
solver_opts = struct('max_iter',1000,'tol',1e-5,'print_level',0);
MPC_problem.solver('ipopt',plugin_opts,solver_opts);
try
    sol = MPC_problem.solve();
    MPC_result = sol.value(MPC_U_input_total);
catch
    MPC_result = MPC_problem.debug.value(MPC_U_input_total);
end
control_applied=MPC_result(1:Nu);
%% 状态更新
arg_state={state_now,control_applied};
state_next_casADi=f_state.call(arg_state);
Power_vector_casADi=f_power.call(arg_state);
state_next=full(state_next_casADi{1});
Power_vector=full(Power_vector_casADi{1});
if state_next(3)>pi
    state_next(3)=state_next(3)-2*pi;
elseif state_next(3)<-pi
    state_next(3)=state_next(3)+2*pi;
end
end
