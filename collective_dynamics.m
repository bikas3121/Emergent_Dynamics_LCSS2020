function [t_col,x_col,x_iso]= collective_dynamics(dwell_time, x_0_in, a_rand, A_cl_1, A_cl_2, t_N)

t_col = [];
x_col = [];
x_iso = [];
 t_span = [0 dwell_time];
x_update = x_0_in;
for i = 1:t_N
    if mod(a_rand(i),2) ==1
        [t1,x1] = ode23s(@(t1,x1) A_cl_1*x1, t_span, x_update);
        t_col = [t_col; t1];
        x_col = [x_col; x1];
        x_iso = [x_iso; x1(:,1:2)];
        len = length (x_col(:,1));
        x_update = x_col(len,:);
    else
        [t2,x2] = ode23s(@(t2,x2) A_cl_2*x2,t_span,x_update);
        t_col = [t_col; t2];
        x_col = [x_col; x2];
        x_iso = [x_iso; x2(:,1:2)];
        len = length(x_col(:,1));
        x_update = x_col(len,:);
    end
    t_span = t_span+ dwell_time;
end
end