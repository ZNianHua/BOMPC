function [index,point_ref] = calc_target_index(state_now,path_ref,Ld,ind_last)     
N=length(path_ref);
n=100;m=100;
distance=zeros(1,N);
begin_flag=1;
end_flag=N;
if nargin<4
    begin_flag=2;
    end_flag=1;
else
    if ind_last-n>0&&ind_last+m<N
        begin_flag=ind_last-n;
        end_flag=ind_last+m-1;
    elseif ind_last-n<1&&ind_last+m<N
        begin_flag=1;
        end_flag=ind_last+m-1;
    elseif ind_last-n>0&&ind_last+m>N
        begin_flag=ind_last-n;
        end_flag=N;
    end
end
if begin_flag>end_flag
    begin_flag=1;
    end_flag=N;
end
for i=begin_flag:end_flag
    distance(i)=norm(state_now(1:2)-path_ref(1:2,i));
end
[~,index]=min(distance(begin_flag:end_flag));
index=index+begin_flag-1;
L=0;
if Ld>=0
    while L<Ld&&index<N
    %     L=L+sqrt(1+tan(path_ref(3,ind))^2)*(path_ref(1,ind+1)-path_ref(1,ind));
        L=L+norm(path_ref(1:2,index+1)-path_ref(1:2,index));
        index=index+1;
    end
else
    index=N;
    while L>Ld&&index>1
    %     L=L+sqrt(1+tan(path_ref(3,ind))^2)*(path_ref(1,ind+1)-path_ref(1,ind));
        L=L-norm(path_ref(1:2,index)-path_ref(1:2,index-1));
        index=index-1;
    end
    if index==1
        index=floor(N/2);
    end
end
point_ref=path_ref(:,index);
    
