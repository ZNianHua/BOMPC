clc;clear;close all;
%% MPC编队测试
%% 参考路径
k_a=20;k_b=10;
path_length=5000;
cx=1/path_length:1000/path_length:1000;
cy=-sin(cx/k_a).*cx/k_b;
cx_last=[0,cx];
cy_last=[0,cy];
delta_x=cx(1:path_length)-cx_last(1:path_length);
delta_y=cy(1:path_length)-cy_last(1:path_length);
ctheta=atan2(delta_y,delta_x);
path_ref=[cx;cy;ctheta];
%% 仿真初始化
T_max=4000;Tstep=0.01;
Nx_4states=4;
Nx_7states=7;
Nu=2;ind_leader=1;
N=3;%跟随者数量
d0=10;%恒间距策略距离
T_follower=200;%跟随者启动时间
phi=0;%路面坡度
uncertain_flag=0.5*rand(7,1);%噪声指示
%% MPC 初始化
Np=30;%预测域
Nc=5;%控制域
MPC_result=zeros(Nu*Nc,N);
Kd_X=1000;
Kd_Y=1000;
Ktheta=10000;
Kv_x=500;
Kv_y=500;
Komega=2500;
Ka=250;
K_Tq=0.1;
K_delta=5000;
K_energy=5e-8+5e-9;
weight_vector=[Kd_X;Kd_Y;Ktheta;Kv_x;Kv_y;Komega;Ka;K_Tq;K_delta;K_energy];
%% 领航者
state_leader=zeros(Nx_7states,T_max+1);
U_leader=zeros(Nu,T_max);
r_tire=0.3;T_power_max=150;i0_all=4.1;eta_all=0.85;
M=1550;a=1.5;b=1.3;Ce=0.01;f=0.015;g=9.8;
Iz=1750;Cx=0.01;Cy=0.06;Cf=-90000;Cr=-70000;hg=0.6;tao=0.11;miu_f=1;miu_r=0;
vehicle_param_4states=[M;a;b;Ce;f;g];
vehicle_param_7states=[M;Iz;a;b;Cx;Cy;Cf;Cr;hg;f;tao;miu_f;miu_r;r_tire;i0_all;eta_all;g];
k=0.5;
target_v=20*ones(1,T_max);%目标车速
%% 跟随者随机初状态
kappa=zeros(T_max,N);
error_state=4;
tracking_error=zeros(error_state,T_max,N);
platoon_error=zeros(T_max,N);
state_follower=zeros(Nx_7states,T_max+1,N);
Power_follower=zeros(4,T_max,N);
ind_platoon=2*ones(N);
U_follower=zeros(Nu,T_max+1,N);
distance=zeros(T_max,N);
vehicle_param_follower=zeros(17,N);
for ind=1:N
    state_follower(1,T_follower,ind)=rand(1)*20-20;
    state_follower(2,T_follower,ind)=rand(1)*40-20;
    state_follower(3,T_follower,ind)=rand(1)*0.6-0.3;
    state_follower(4,T_follower,ind)=rand(1)*10;
    coi=0.8+0.5*rand(1);
    vehicle_param_follower(1,ind)=M*coi;
    vehicle_param_follower(2,ind)=Iz*coi;
    vehicle_param_follower(3,ind)=a*coi;
    vehicle_param_follower(4,ind)=b*coi;
    vehicle_param_follower(5,ind)=Cx*coi;
    vehicle_param_follower(6,ind)=Cy*coi;
    vehicle_param_follower(7,ind)=Cf*coi;
    vehicle_param_follower(8,ind)=Cr*coi;
    vehicle_param_follower(9,ind)=hg*coi;
    vehicle_param_follower(10,ind)=f*coi;
    vehicle_param_follower(11,ind)=tao*coi;
    vehicle_param_follower(12,ind)=miu_f*coi;
    vehicle_param_follower(13,ind)=1-vehicle_param_follower(12,ind);
    vehicle_param_follower(14,ind)=r_tire*coi;
    vehicle_param_follower(15,ind)=i0_all;
    vehicle_param_follower(16,ind)=eta_all;
    vehicle_param_follower(17,ind)=g;
end
%% 仿真过程
tic
for i=1:T_max
    %% 领航者控制：纯跟踪控制+P控制
    Ld=k*state_leader(4,i)-1;
    [ind_leader,point_ref] = calc_target_index(state_leader(:,i),path_ref,Ld,ind_leader);%找最近点
    [delta] = pure_pursuit_control(state_leader(:,i),point_ref,a+b,Ld);%纯跟踪
    Tq_leader=50*(target_v(i)-state_leader(4,i));
    U_leader(:,i)=[Tq_leader;delta];
    state_leader(:,i+1)= VehicleModel(state_leader(:,i),vehicle_param_7states,U_leader(:,i),phi,Tstep);
    %% 跟随者：侧偏力模型，MPC耦合控制
    %% 跟随者时间启动
    if i>T_follower-1
        for j=1:N
            platoon_d=-j*d0;
            %% 跟随者控制：
            vehicle_param_follower(:,j)=vehicle_param_7states;
            sequence=state_leader(:,1:i);
%             if j==1
%                 sequence=state_leader(:,1:i);
%             else
%                 sequence=state_follower(:,1:i,j-1);
%             end
            K=[3 3];
            M_follower=vehicle_param_follower(1,j);
            Ca_follower=vehicle_param_follower(5,j);
            f_follower=vehicle_param_follower(10,j);
            [a_i,omega_i,kappa(i,j)] = controller_compare(state_follower(:,i,j),sequence,j*d0,Tstep,kappa(i-1,j),K);
            [T_i,delta_i] = comparedInputConvert(a_i,omega_i,vehicle_param_follower(:,j),state_follower(:,i,j),phi);
            U_follower(:,i,j)=[T_i,delta_i];
            [state_follower(:,i+1,j),Power_follower(:,i,j)]=VehicleModel(state_follower(:,i,j),vehicle_param_follower(:,j),U_follower(:,i,j),phi,Tstep);
            [distance(i,j)] = calculate_S(state_leader(:,1:i),state_follower(:,i,j));%计算曲线坐标
            platoon_error(i,j)=distance(i,j)-j*d0;
%             if rem(i,10)==1
%             scatter(state_leader(1,i),state_leader(2,i),'+')
%             hold on
%             scatter(state_follower(1,i+1,j),state_follower(2,i+1,j));
%             pause(0.1)
%             end 
        end
    end
end
toc
%% 绘图
subplot(4,2,[1,3])
plot(state_leader(1,:),state_leader(2,:))
for ii=1:N
    hold on
    plot(state_follower(1,T_follower:end,ii),state_follower(2,T_follower:end,ii))
end
legend('leader','follower_1','follower_2','follower_3','follower_4')
subplot(4,2,[2,4])
plot(state_leader(4,:));
for iii=1:N
    hold on
    plot(state_follower(4,:,iii))
end
legend('leader','follower_1','follower_2','follower_3','follower_4')
subplot(4,2,[5,7])
for iiii=1:N
    plot(distance(:,iiii))
    hold on
end
legend('follower_1','follower_2','follower_3','follower_4')
subplot(4,2,6)
for iiii=1:N
    plot(abs(sum(Power_follower(:,100:end,iiii),1)));
    hold on
end
legend('follower_1','follower_2','follower_3','follower_4')
subplot(4,2,8)
Power_v=Energy_XandYandYAW(state_follower,U_follower,T_max,T_follower,N,vehicle_param_follower);
for iiii=1:N
    plot(Power_v(:,iiii))
    hold on
end
legend('follower_1','follower_2','follower_3','follower_4')
