function [t_mf,x_mf]= mean_field_dynamics(dwell_time, x_0_mf, t_N, a_rand, mf_1,mf_2, J_1, J_2)
x_mf =[];
t_mf= [];
 t_span = [0 dwell_time];
x_mf_up = x_0_mf;
for i = 1:t_N
    if mod(a_rand(i),2) ==1
        [t5,x5] = ode23s(@(t5,x5) mf_1*x5, t_span, x_mf_up);
        t_mf = [t_mf; t5];
        x_mf = [x_mf; x5];
        len = length (x_mf(:,1));
        x_mf_up = (J_1* x_mf(len,:)')';
    else
        [t6,x6] = ode23s(@(t6,x6) mf_2*x6,t_span,x_mf_up);
        t_mf = [t_mf; t6];
        x_mf = [x_mf; x6];
        len = length(x_mf(:,1));
        x_mf_up = (J_2*x_mf(len,:)')';
    end
    t_span = t_span+ dwell_time;
end
