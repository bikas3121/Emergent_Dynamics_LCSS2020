%% 
clc
clear 
close all
N= 30; %number of agents

nx = 2;

%% Input
G =100;
dwell_time = 0.2;
preferred_no_of_unstable_modes = 500 ;
preferred_no_of_stable_modes =  0.45 * preferred_no_of_unstable_modes ;
total_switchings = preferred_no_of_stable_modes + preferred_no_of_unstable_modes;
t_N = total_switchings;

%% State Matrices
load data/A4.mat

%% Graph 
format short
load data/s1_30.mat
load data/t1_30.mat
load data/s2_30.mat
load data/t2_30.mat

%% Laplacian Generation
L1 = lap_gen_indeg(N,s1_30,t1_30);
L2 = lap_gen_indeg(N,s2_30,t2_30);

%% Eigenvvectors of the laplacians
[V1, D1, v1] = normalize_eigenvector(L1);
[V2, D2, v2] = normalize_eigenvector(L2);

%Removing the zero eigenvalue from the D1
for i = 1:N
    if D1(i,i) == 0
        p =i;
    end
end
D1(p,:) = [];
D1(:,p) = [];

for i = 1:N
    if D2(i,i) == 0
        p2 =i;
    end
end
D2(p2,:) = [];
D2(:,p2) = [];

%% initial conditions
%x_0_in = -1 + (1+1)*rand(2*N,1);% Initial Values
load data/x_0_in.mat;

%% Collective Dynamics
A_cl_1  =  (A - G*kron(L1,eye(nx))); 
A_cl_2  =  (A - G*kron(L2,eye(nx)));

%% MeanField and Error Dynamics
A_0_1 = kron(v1,eye(nx))*A*kron(ones(N,1),eye(nx));
A_0_2 = kron(v2,eye(nx))*A*kron(ones(N,1),eye(nx));
B_1_1 = kron(v1,eye(nx))*A*kron(V1,eye(nx));
B_1_2 = kron(v2,eye(nx))*A*kron(V2,eye(nx));
B_2_1 = kron(pinv(V1),eye(nx))*A*kron(ones(N,1),eye(nx));
B_2_2 = kron(pinv(V2),eye(nx))*A*kron(ones(N,1),eye(nx));
B_3_1 = kron(pinv(V1),eye(nx))*A*kron(V1,eye(nx));
B_3_2 = kron(pinv(V2),eye(nx))*A*kron(V2,eye(nx));
A_22_1 = -kron(D1,eye(nx));
A_22_2 = -kron(D2,eye(nx));
mf_1 = [A_0_1 B_1_1; B_2_1 (G*A_22_1+B_3_1)];
mf_2 = [A_0_2 B_1_2; B_2_2 (G*A_22_2+B_3_2)];


eig1 = eig(A_0_1+ A_0_1');
eig2 = eig(A_0_2+A_0_2');
ratio = eig2/eig1;


%% Jump Dynamics
% Jump from graph 1 to graph 2
J_11_1 = eye(2);
J_12_1 = kron(v2*V1,eye(2)) ;
J_21_1 = zeros(2*N-2,2);
J_22_1 = kron(pinv(V2)*V1,eye(2));
J_1 = [J_11_1 J_12_1; J_21_1 J_22_1];

% Jump from graph 2 to graph 1
J_11_2 = eye(2);
J_12_2 = kron(v1*V2,eye(2)) ;
J_21_2 = zeros(2*N-2,2);
J_22_2 = kron(pinv(V1)*V2,eye(2));
J_2 = [J_11_2 J_12_2; J_21_2 J_22_2];

%% Initial values for Emergent Dynamics and Error over G1
x_0_avg1 = kron(v1,eye(nx))*x_0_in;
%e_0_1 = x_0_in - kron(ones(N,1),x_0_avg1);
e_v_1 = kron(pinv(V1),eye(nx))*x_0_in;
x_0_mf1 = [x_0_avg1; e_v_1];

%% Initial values for Emergent Dynamics and Error over G2
x_0_avg2 = kron(v2,eye(nx))*x_0_in;
%e_0_2 = x_0_in - kron(ones(N,1),x_0_avg2);
e_v_2 = kron(pinv(V2),eye(nx))*x_0_in;
x_0_mf2 = [x_0_avg2; e_v_2];

%%  Random Switching sequence
%3940/3000


%% Random Switching Pattern
a1 = [];
for i = 1: preferred_no_of_stable_modes
    a1 = [a1 1];
end

a2 = [];
for i = 1: preferred_no_of_unstable_modes
    a2 = [a2 2];
end
a = [a1 a2];
 a_rand = a(randperm(length(a)));
%% Fixed Switching Pattern

% load data/a_rand.mat
% l= length(a_rand);
% if l>= t_N
%     a_rand = a_rand;
% else
%     dif = t_N-l;
%     dif2 = ceil(dif/l);
%     for i = 1  : dif2
%         a_rand = [a_rand a_rand];
%     end
% end

%% Initial Value Allocation
if a_rand(1) == 1
    x_0_avg = x_0_avg1;
    x_0_mf = x_0_mf1;
else
    x_0_avg = x_0_avg2;
    x_0_mf = x_0_mf2;
end

%% Switching Collective Dynamics
t_col = [];
x_col = [];
 t_span = [0 dwell_time];
x_update = x_0_in;
for i = 1:t_N
    if mod(a_rand(i),2) ==1
        [t1,x1] = ode23s(@(t1,x1) A_cl_1*x1, t_span, x_update);
        t_col = [t_col; t1];
        x_col = [x_col; x1];
        len = length (x_col(:,1));
        x_update = x_col(len,:);
    else
        [t2,x2] = ode23s(@(t2,x2) A_cl_2*x2,t_span,x_update);
        t_col = [t_col; t2];
        x_col = [x_col; x2];
        len = length(x_col(:,1));
        x_update = x_col(len,:);
    end
    t_span = t_span+ dwell_time;
end


len_col = length(t_col);
x_col1 = [];
x_col2 = [];

for i = 1:nx*N
    if mod(i,2)==1 
        x_col1 = [x_col1 x_col(:,i)];
    else
        x_col2 = [x_col2 x_col(:,i)];
    end
end




%% Emergent Dynamics
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


%% Mean-Field Dynamics
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

%% Plots
figure
subplot(2,1,1);
 plot(t_col,x_col(:,1))
 hold on 
  plot(t_em,x_em(:,1),'--', 'color','green')
  hold on
 plot(t_mf,x_mf(:,1),'-.','color', 'red')
xlabel('Time')
legend('x(1)','x_{e}(1)','x_{s}(1)')
grid on

subplot(2,1,2); 
 plot(t_col,x_col(:,2))
 hold on
  plot(t_em,x_em(:,2),'--', 'color','green')
  hold on
plot(t_mf,x_mf(:,2),'-.','color', 'red')
legend('x(2)','x_{e}(2)','x_{s}(2)')
xlabel('Time')
grid on

