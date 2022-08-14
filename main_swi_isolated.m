%% 
clc
clear 
close all
N = 30; %number of agents
nx = 2; %dimension of the system
G =100; % interconnection strength,gamma.
dwell_time = 10; %dwell time in each mode.

%% Input
% we choose the following numerical value for the simulation presented in
% the paper. 



%% State Matrices
load data_iso/A.mat %state matrices of the agents.

%% Graph 
% Generation of the strongly connected directed graph. 
% File s1_30.mat and t1_30.mat contains the edge sets for graph G1.
% File s2_30.mat and t2_30.mat contains the edge sets for graph G2.
format short
load data_iso/s1_30.mat
load data_iso/t1_30.mat
load data_iso/s2_30.mat
load data_iso/t2_30.mat

G1 = digraph(s1_30, t1_30); %generates the graph G1
G2 = digraph(s2_30, t2_30); %generates the graph G2
% figure 
% plot(G1) %plots the graph G1.
% 
% figure 
% plot(G2) %plots the graph G2.

%% Laplacian Generation
% The function 'lap_gen_indeg' generate the Laplacian of the graph where
% number of nodes(N), head of the edge (s) and tail of the edge in the
% graph (t) are given.
L1 = lap_gen_indeg(N,s1_30,t1_30); %Laplacian of the graph G1 
L2 = lap_gen_indeg(N,s2_30,t2_30); %Laplacian of the graph G2

%% Eigenvector of the Laplacians with isolated node
[V1,D1,W1] = eig(L1); %Calculates the eigenvalue and eigenvectors of Laplacian L1
[V2,D2,W2] = eig(L2); %Calculates the eigenvalue and eigenvectors of Laplacian L2
% This takes only the real  part of the eigenvalue and eigenvectors.
V1 = real(V1);
D1 = real(D1);
W1 = real(W1);
V2 = real(V2);
D2 = real(D2);
W2 = real(W2);

% making the right eigenvector corresponding to the zero eigenvalue as
% vector of ones. 
V1 = V1/0.1857;
V2 = V2/0.1857;
% normalizing the left eigenvectors of the laplacian corresponding to the
% zero eigenvalue.
sm1 = sum(W1(:,1));
W1 = W1/sm1;
sm2 = sum(W2(:,3));
W2 = W2/sm2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%left-eigenvector corresponding to the zero eigenvalue
v1 = W1(:,1)'; %left-eigenvector corresponding to the zero eigenvalue of L1
v2 = W2(:,3)'; %left-eigenvector corresponding to the zero eigenvalue of L2

%other generalized eigenvectors
V1(:,1) = []; %right eigenvector of the L1 except the vector of ones.
V2(:,3) = []; %right eigenvector of the L2 except the vector of ones.

%Removing the zero eigenvalue from the D1
% since the eigenvalues that should be zero are not exactly zero in MATLAB,
% hence we replace very small number with the exact zero. 
% In the diagonal matrix with the eigenvalue of the Laplacian L1, the first
% eigenvalue is zero.
for i = 1:N
    if D1(i,i) <= 0.00000000001
        D1(i,i) =0;
    end
end
% remove the row and column corresponding the zero eigenvalue.
D1(1,:) = [];
D1(:,1) = [];
% In the diagonal matrix with the eigenvalue of the Laplacian L2, the third
% eigenvalue is zero.
for i = 1:N
    if D2(i,i) <= 0.0000000001
        D2(i,i) = 0 ;
    end
end
% remove the row and column corresponding the zero eigenvalue.
D2(3,:) = [];
D2(:,3) = [];

%% initial conditions
%x_0_in = -100 + (100+100)*rand(2*N,1);% this generates the random initial
%values.

load data_iso/x_0_in.mat; % this initial values were used in the simulation result presented in the paper. 

%% Collective Dynamics
% The calculates the state matrices for the collective closed-loop
% dynamics, equation 4 in the paper. 
A_cl_1  =  (A - G*kron(L1,eye(nx))); 
A_cl_2  =  (A - G*kron(L2,eye(nx)));

%% MeanField and Error Dynamics
% The following matrices are the mean-field and error dynamics with graph
% G1
%refer to equation 20 in the paper for the graph G1
A_0_1 = kron(v1,eye(nx))*A*kron(ones(N,1),eye(nx));
B_1_1 = kron(v1,eye(nx))*A*kron(V1,eye(nx));
B_2_1 = kron(pinv(V1),eye(nx))*A*kron(ones(N,1),eye(nx));
B_3_1 = kron(pinv(V1),eye(nx))*A*kron(V1,eye(nx));
A_22_1 = -kron(D1,eye(nx));
mf_1 = [A_0_1 B_1_1; B_2_1 (G*A_22_1+B_3_1)];

% The following matrices are the mean-field and error dynamics with graph
% G1
%refer to equation 20 in the paper for the graph G2
A_0_2 = kron(v2,eye(nx))*A*kron(ones(N,1),eye(nx));
B_1_2 = kron(v2,eye(nx))*A*kron(V2,eye(nx));
B_2_2 = kron(pinv(V2),eye(nx))*A*kron(ones(N,1),eye(nx));
B_3_2 = kron(pinv(V2),eye(nx))*A*kron(V2,eye(nx));
A_22_2 = -kron(D2,eye(nx));
mf_2 = [A_0_2 B_1_2; B_2_2 (G*A_22_2+B_3_2)];


%% Jump Dynamics

% Jump from graph 1 to graph 2
%refer to equation 22 in the paper. 
J_11_1 = eye(2);
J_12_1 = kron(v2*V1,eye(2)) ;
J_21_1 = zeros(2*N-2,2);
J_22_1 = kron(pinv(V2)*V1,eye(2));
J_1 = [J_11_1 J_12_1; J_21_1 J_22_1];

% Jump from graph 2 to graph 1
%refer to equation 22 in the paper. 
J_11_2 = eye(2);
J_12_2 = kron(v1*V2,eye(2)) ;
J_21_2 = zeros(2*N-2,2);
J_22_2 = kron(pinv(V1)*V2,eye(2));
J_2 = [J_11_2 J_12_2; J_21_2 J_22_2];

%% Calulate the eigenvalue of the emergent dynamics for the verifying the boundedness condition. 
% refer to equation 27 in the paper.
eig1 = 0.5*max(eig(A_0_1+ A_0_1')); %for graph G1
eig2 = 0.5*max(eig(A_0_2+A_0_2')); %for graph G2

%% Boundedness condition
% to satisfy the boundedness condition, i.e., the Lemma 2, the following
% ratio should be satisfied. 
% refer to the equation 29 in the paper, with the dwell time =10 sec(in the
% simulation in the paper).
ratio = eig2/eig1;
% this show that if we choose 10 unstable nodes (with dwell time =10)  then
% the number of  stable nodes should be ratio*(no of unstable nodes)
preferred_no_of_unstable_modes = 10 ; 
preferred_no_of_stable_modes =  round(0.8779 * preferred_no_of_unstable_modes) ;
total_switchings = preferred_no_of_stable_modes + preferred_no_of_unstable_modes;
t_N = total_switchings; %total switching is the number of the stable and unstable nodes. 


%% Initial values for Emergent Dynamics and Error 
% due to the coordinate transformation the initial value use for the
% collective dynamics should be transformed using the same transformation
% used for transforming the system dynamics into mean-field and error
% dynamics. 
 %see the equations 10 and 12 from the paper. 
%For graph G1
x_0_avg1 = kron(v1,eye(nx))*x_0_in;
e_v_1 = kron(pinv(V1),eye(nx))*x_0_in;
x_0_mf1 = [x_0_avg1; e_v_1];

% For graph  G2
x_0_avg2 = kron(v2,eye(nx))*x_0_in;
e_v_2 = kron(pinv(V2),eye(nx))*x_0_in;
x_0_mf2 = [x_0_avg2; e_v_2];



%% Random Switching Pattern
% This generates the random switching pattern. 
a1 = [];
for i = 1: preferred_no_of_stable_modes
    a1 = [a1 1];
end

a2 = [];
for i = 1: preferred_no_of_unstable_modes
    a2 = [a2 2];
end
a = [a1 a2];
%a_rand = a(randperm(length(a)));
load data_iso/a_rand1.mat % this switching pattern in used for the simulation in the paper. 


%% Initial Value Allocation
% this choose the initial value of based in switching generated by the
% 'Random Switching Pattern'.
% the initial value for the transformed system depends on the graph. 
if a_rand(1) == 1
    x_0_avg = x_0_avg1;
    x_0_mf = x_0_mf1;
else
    x_0_avg = x_0_avg2;
    x_0_mf = x_0_mf2;
end

%% Switching Collective Dynamics

[t_col,x_col,x_iso]= collective_dynamics(dwell_time, x_0_in, a_rand, A_cl_1, A_cl_2, t_N);

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
[t_em,x_em]= emergent_dynamics(dwell_time, x_0_avg, t_N, a_rand, A_0_1, A_0_2);

%% Mean-Field Dynamics
[t_mf,x_mf]= mean_field_dynamics(dwell_time, x_0_mf, t_N, a_rand, mf_1,mf_2, J_1, J_2);


%%
p(1) = plot(t_col,x_iso(:,1));
hold on 
 plot(t_col,x_col1(:,2:30));
 hold on
 p(2) = plot(t_col,x_col1(:,2));
hold on 
p(3) = plot(t_em,x_em(:,1),'--', 'color','green');
hold on
p(4) = plot(t_mf,x_mf(:,1),'-.','color', 'red');

hold on

% display the switching sequence
dec = 10*rem(t_N,1);
if dec >=5 
    t_N = ceil(t_N);
else 
    t_N = floor(t_N);
end
swi = [];
for i = 1:t_N
    a =  i*dwell_time;
    swi = [swi a];
end
for i = 1:length(swi)
    xline(swi(i))
    hold on 
end
hold on 
xlabel('Time (sec)')
legend(p(1:4),'x_{isolated}(1)', 'x(1)', 'x_e (1)', 'x_{s}(1)')
grid on

% boxed plot 
axes('position',[.65 .175 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (t_col < 0.006) & (t_col >=0); % range of t near perturbation
plot(t_col(indexOfInterest),x_col1(indexOfInterest,:)); % plot on new axes
xticks([0 0.002 0.004 0.006])
xticklabels({'0','0.002','0.004','0.006'})
%axis fill
grid on
ylim([-50 50])





