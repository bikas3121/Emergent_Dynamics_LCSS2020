function [t_em,x_em]= emergent_dynamics(dwell_time, x_0_avg, t_N, a_rand, A_0_1, A_0_2)
% This function solves the emergent dynamics, i.e. the reduced order slow
% dynamics obtained after the time-scale separation
% INPUTS:
  % dwell_time, 
  % x_0_avg -  initial value to the mean-field dynamics obtained after the
  % transfromation
  % t_N -  Simulation time
  % a_rand - switching mode
  % A_0_1 and A_0_2 -  Closed loop network dynamics due to Graph G1 and G2,
  % respectively
 
x_em =[];
t_em= [];
 t_span = [0 dwell_time];
x_em_up = x_0_avg;
for i = 1:t_N
    if mod(a_rand(i),2) ==1
        [t3,x3] = ode23s(@(t3,x3) A_0_1*x3, t_span, x_em_up);
        t_em = [t_em; t3];
        x_em = [x_em; x3];
        len = length (x_em(:,1));
        x_em_up = x_em(len,:);
    else
        [t4,x4] = ode23s(@(t4,x4) A_0_2*x4,t_span,x_em_up);
        t_em = [t_em; t4];
        x_em = [x_em; x4];
        len = length(x_em(:,1));
        x_em_up = x_em(len,:);
    end
    t_span = t_span + dwell_time;
end
end
