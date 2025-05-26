function [L] = calculate_S(state_leader,state_now)
[~,m]=size(state_leader);
L=0;
while state_leader(1,m)>state_now(1)&&m>1
    L=L+norm(state_leader(1:2,m)-state_leader(1:2,m-1));
    m=m-1;
end
end

