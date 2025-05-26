function [reference_path,original_path,N_frequency] = path_generate(all_path,target_ind,Tstep,Np)
rest_path=all_path(:,target_ind:end);
[m,~]=size(rest_path);
path_out=zeros(m,Np);
path_out(:,1)=rest_path(:,1);
%% 信息频率观测
delta_x=(rest_path(1,2:end)-rest_path(1,1:end-1));
delta_x=[delta_x,delta_x(end)];
Tstep_series=delta_x./(rest_path(4,:).*cos(rest_path(3,:)));
% Tstep_accept=mean(Tstep_series);
Tstep_accept=mode(Tstep_series(Tstep_series>0));
%% 步长对比
% N_frequency=round(Tstep_accept/Tstep);
N_frequency=1;
%% 均匀变化值解算
dstate=(rest_path(:,2:N_frequency:end)-rest_path(:,1:N_frequency:end-1))./Tstep_accept;
[~,n]=size(dstate);
% dstate=[dstate,dstate(:,end)];
%% 线性插值
j=1;
delta_state=dstate(:,j);
for i=1:Np
    if i*Tstep>=j*Tstep_accept
        j=j+1;
        if j>n-1
            j=n-1;
        end
        delta_state=dstate(:,j);
    end
    path_out(:,i+1)=path_out(:,i)+delta_state*Tstep;
end
original_path=rest_path(:,1:j);
reference_path=path_out(:,2:end);
end