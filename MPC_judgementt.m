function J_cost = MPC_judgementt(obj_param,state_leader,platoon_param,Time_param,User_param)
%% 数据提取
fix_param=platoon_param{1};
Nx=fix_param(1);
Nu=fix_param(2);
N=fix_param(3);
d0=fix_param(4);%恒间距策略距离
vehicle_param_follower=platoon_param{2};

T_max=Time_param(1);
T_follower=Time_param(2);%跟随者启动时间
Tstep=Time_param(3);

phi=0;%路面坡度
%% MPC 初始化
% if length(obj_param)<5
%     weight_vector=obj_param{1};
%     MPC_hori=obj_param{2};
%     Np=MPC_hori(1);
%     Nc=MPC_hori(1);
% else
%     Np=30;%预测域
%     Nc=5;%控制域
%     weight_vector=obj_param;
% end

% Kd_X=obj_param.Kd_X;
% Kd_Y=1000;
% Ktheta=10000;
% Kv_x=300;
% Kv_y=300;
% Komega=2000;
% Ka=250;
% K_Tq=0.1;
% K_delta=10000;
[m,n]=size(obj_param);
if max([m,n])<31
    Np=30;%预测域
    Nc=5;%控制域
    for j=1:N
        switch j
            case 1
                weight_vector_1=[obj_param.Kd_X_1;obj_param.Kd_Y_1;obj_param.Ktheta_1;
                    obj_param.Kv_x_1;obj_param.Kv_y_1;obj_param.Komega_1;
                    obj_param.Ka_1;obj_param.K_Tq_1;obj_param.K_delta_1;obj_param.K_energy_1];
            case 2
                weight_vector_2=[obj_param.Kd_X_2;obj_param.Kd_Y_2;obj_param.Ktheta_2;
                    obj_param.Kv_x_2;obj_param.Kv_y_2;obj_param.Komega_2;
                    obj_param.Ka_2;obj_param.K_Tq_2;obj_param.K_delta_2;obj_param.K_energy_1];
            case 3
                weight_vector_3=[obj_param.Kd_X_3;obj_param.Kd_Y_3;obj_param.Ktheta_3;
                    obj_param.Kv_x_3;obj_param.Kv_y_1;obj_param.Komega_3;
                    obj_param.Ka_3;obj_param.K_Tq_3;obj_param.K_delta_3;obj_param.K_energy_1];
        end
    end
    switch N
        case 1
            weight_vector=weight_vector_1;
        case 2
            weight_vector=[weight_vector_1,weight_vector_2];
        case 3
            weight_vector=[weight_vector_1,weight_vector_2,weight_vector_3];
    end
else
    Np=obj_param.Np;%预测域
    Nc=obj_param.Nc;%控制域
    for j=1:N
        switch j
            case 1
                weight_vector_1=[obj_param.Kd_X_1;obj_param.Kd_Y_1;obj_param.Ktheta_1;
                    obj_param.Kv_x_1;obj_param.Kv_y_1;obj_param.Komega_1;
                    obj_param.Ka_1;obj_param.K_Tq_1;obj_param.K_delta_1;obj_param.K_energy_1];
            case 2
                weight_vector_2=[obj_param.Kd_X_2;obj_param.Kd_Y_2;obj_param.Ktheta_2;
                    obj_param.Kv_x_2;obj_param.Kv_y_2;obj_param.Komega_2;
                    obj_param.Ka_2;obj_param.K_Tq_2;obj_param.K_delta_2;obj_param.K_energy_1];
            case 3
                weight_vector_3=[obj_param.Kd_X_3;obj_param.Kd_Y_3;obj_param.Ktheta_3;
                    obj_param.Kv_x_3;obj_param.Kv_y_1;obj_param.Komega_3;
                    obj_param.Ka_3;obj_param.K_Tq_3;obj_param.K_delta_3;obj_param.K_energy_1];
        end
    end
    switch N
        case 1
            weight_vector=weight_vector_1;
        case 2
            weight_vector=[weight_vector_1,weight_vector_2];
        case 3
            weight_vector=[weight_vector_1,weight_vector_2,weight_vector_3];
    end
end
MPC_result=zeros(Nu*Nc,N);
%% 跟随者随机初状态
error_state=User_param(1);%指定计入误差判定的前n个状态
error_judge_matrix=diag(User_param(2:1+error_state).^2);%指定前n个状态的权值（未统一数量级）
platoon_judge_factor=User_param(2+error_state).^2;%指定编队间距误差的权值（未统一数量级）
power_judge_factor=User_param(3+error_state).^2;%指定功率权值（未统一数量级）

tracking_error=zeros(T_max,N);
platoon_error_single=zeros(T_max,N);
state_follower=zeros(Nx,T_max+1,N);
power_follower=zeros(T_max,N);
ind_platoon=2*ones(N);
U_follower=zeros(Nu,T_max+1,N);
distance=zeros(T_max,N);
for ind=1:N
    state_follower(1,T_follower,ind)=rand(1)*20-20;
    state_follower(2,T_follower,ind)=rand(1)*40-20;
    state_follower(3,T_follower,ind)=rand(1)*0.6-0.3;
    state_follower(4,T_follower,ind)=rand(1)*10;
end
%% 仿真过程
record_judge_flag=ones(N,1);
recorded_stable_time=T_follower.*ones(N,1);
ind_control_error=ones(N,1);
recorded_stable_time_control_error=ones(N,1);
cost_total_single_vehicle=zeros(N,1);
for i=1:T_max
    %% 跟随者：侧偏力模型，MPC耦合控制
    %% 跟随者时间启动
    if i>T_follower-1
        for j=1:N
            %% 跟随者控制：
            if rem(i,100)==1%减低控制频率
                [ind_platoon(j),~] = calc_target_index(state_follower(:,i,j),state_leader(:,1:i),-j*d0,ind_platoon(j));%找对应跟随者的期望位置
                [path_out,~,~] = path_generate(state_leader(:,1:i),ind_platoon(j),Tstep,Np);%线性插值，不同频率下生成跟踪轨迹
                [state_follower(:,i+1,j),U_follower(:,i,j),MPC_result(:,j),power_vector_follower] = MPCPlatoonAndUpdate(state_follower(:,i,j),path_out,vehicle_param_follower(:,j),phi,Tstep,[Np;Nc],weight_vector(:,j),MPC_result(:,j));   
                tracking_error(ind_control_error(j),j)=(path_out(1:error_state,1)-state_follower(1:error_state,i+1,j))'*error_judge_matrix*(path_out(1:error_state,1)-state_follower(1:error_state,i+1,j));
                ind_control_error(j)=ind_control_error(j)+1;
            else
                U_follower(:,i,j)= U_follower(:,i-1,j);
                [state_follower(:,i+1,j),power_vector_follower]=VehicleModel(state_follower(:,i,j),vehicle_param_follower(:,j),U_follower(:,i,j),phi,Tstep);
            end
            [distance(i,j)] = calculate_S(state_leader(:,1:i),state_follower(:,i,j));%计算曲线坐标
            platoon_error_single(i,j)=distance(i,j)-j*d0;
            power_follower(i,j)=power_vector_follower'*eye(4)*power_vector_follower;
            if ind_control_error(j)>1 && tracking_error(ind_control_error(j)-1,j)>0 && tracking_error(ind_control_error(j)-1,j)<15 && record_judge_flag(j)==1
                recorded_stable_time(j)=i;
                recorded_stable_time_control_error(j)=ind_control_error(j);
                record_judge_flag(j)=0;
            end
        end
    end
end
%% BO代价函数计算
for j=1:N
    loss_power=power_judge_factor*power_follower(recorded_stable_time(j):T_max,j)/(T_max-recorded_stable_time(j)+1);
    control_error=tracking_error(recorded_stable_time_control_error(j):ind_control_error(j)-1,j)/(ind_control_error(j)-recorded_stable_time_control_error(j));%MPC控制误差
    platoon_error=platoon_judge_factor*platoon_error_single(recorded_stable_time(j):T_max,j).^2/(T_max-recorded_stable_time(j)+1);%编队间距误差
    cost_total_single_vehicle(j)=sum(loss_power)+sum(control_error)+sum(platoon_error);
end
J_cost=sum(cost_total_single_vehicle);





