clc;clear;close all;
%% 数据初始化
N_repeat = 5;
N_iter=10;
Nx=7;
Nu=2;
N=3;
d0=10;
T_max=4000;
T_follower=20;
Tstep=0.01;
fix_param=[Nx;Nu;N;d0];
r_tire=0.3;T_power_max=150;i0_all=4.1;eta_all=0.85;
M=1550;a=1.5;b=1.3;Ce=0.01;f=0.015;g=9.8;
Iz=1750;Cx=0.01;Cy=0.06;Cf=-90000;Cr=-70000;hg=0.6;tao=0.08;miu_f=1;miu_r=0;
vehicle_param_follower=zeros(17,N);
for ind=1:N
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
platoon_param={fix_param;vehicle_param_follower};
Time_param=[T_max;T_follower;Tstep];
file_name_part1='leader_traj_';
file_name_part2='.mat';
% state_leader=zeros(Nx,T_max+1);
error_state=4;%指定计入误差判定的前n个状态
kd_x=1;
kd_y=1;
ktheta=10;
kv=0.5;
platoon_judge_factor=1;%指定编队间距误差的权值（未统一数量级）
power_judge_factor=1e-4;%指定功率权值（未统一数量级）
User_param=[error_state;kd_x;kd_y;ktheta;kv;platoon_judge_factor;power_judge_factor];
%% 贝叶斯优化
weight_vector_name={'Kd_X';'Kd_Y';'Ktheta';'Kv_x';'Kv_y';'Komega';'Ka';'K_energy';'K_Tq';'K_delta'};
num_weight=length(weight_vector_name);
for i=1:N
    for j=1:num_weight
        variable_type='integer';
        switch j
            case 1
                lb=500;ub=1500;
            case 2
                lb=500;ub=1500;
            case 3
                lb=5000;ub=15000;
            case 4
                lb=250;ub=750;
            case 5
                lb=250;ub=750;
            case 6
                lb=1250;ub=3750;
            case 7
                lb=125;ub=325;
            case 8
                lb=1e-8;ub=1e-7;
                variable_type='real';
            case 9
                lb=0.05;ub=0.15;
                variable_type='real';
            case 10
                lb=2590;ub=7500;
        end
        xbo_weight(num_weight*i-num_weight+j)=optimizableVariable([weight_vector_name{j},'_',num2str(i)],[lb ub],'Type',variable_type);
    end
end
xbo_Np=optimizableVariable('Np',[10 30],'Type','real');
xbo_Nc=optimizableVariable('Nc',[1 10],'Type','integer');
% xbo=[xbo_weight,xbo_Np,xbo_Nc];
xbo=xbo_weight;

tic
saved_results = cell(N_repeat,1);
for i = 1:N_repeat
    file_name=[file_name_part1,num2str(i),file_name_part2];
    load(file_name);
    fun = @(x) MPC_judgementt(x,state_leader,platoon_param,Time_param,User_param);
    results = bayesopt(fun,xbo,...
        'AcquisitionFunctionName', 'expected-improvement',...
        'IsObjectiveDeterministic', 0,...
        'ExplorationRatio', 0.5,...
        'GPActiveSetSize', 300,...
        'UseParallel', false,...
        'MaxObjectiveEvaluations', N_iter,...
        'PlotFcn', [], ...
        'InitialX', [], ...
        'InitialObjective', [], ...
        'InitialConstraintViolations', [], ...
        'NumSeedPoints', length(xbo)+1);
    saved_results{i} = results;
end
toc
save("BO_result_5.mat")
