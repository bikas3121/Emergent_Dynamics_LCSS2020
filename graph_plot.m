%function [s,t,n_bin,pop_node] = graph_plot (NumberOfNodes)
clc
clear
close all
%% Generate the permutation of the nodes.
perm = []; 
NumberOfNodes = 30;
for i =1:NumberOfNodes
    for j = 1:NumberOfNodes
        perm = [perm;i,j;];
    end
end
 perm(diff(perm,[],2)==0, :) = []; %remove self-loops
 
 %% Generate the edges
 len_perm = length(perm); 
 % To make some node popular we chose the nodes in random and make the
 % number of edges to and from these nodes more than the average. 
 num_pop_node = ceil(0.4*NumberOfNodes);
 pop_node = sort(randperm(NumberOfNodes,num_pop_node),'ascend'); % generate random popular nodes and sort them

 
%% This generates the edges for the graph while considering the popular edges. 
percent_of_edges = floor(0.2*NumberOfNodes);
pos_edges = randi([1,percent_of_edges],NumberOfNodes,1)';      %number of possible edges for each node.
%num_edge_pop_nodes = NumberOfNodes*0.3  ;          % this is the number of edges assigned to the popular node.
num_edge_pop_nodes = randi([ceil(0.2*NumberOfNodes),ceil(NumberOfNodes*0.3)],num_pop_node,1)';
% This assigns the number of edges for the popular nodes. 
for i = 1: NumberOfNodes
    for j = 1:length(pop_node)
             if i == pop_node(j)
                 pos_edges(i) = num_edge_pop_nodes(j);
             end
    end
end


%% Generates the graph edges based on the above construction.
set_edges= [];
for i = 1:length(pos_edges)
   a = pos_edges(i);
   p = (NumberOfNodes-1)*i;
   if i == 1
        pos = randperm(NumberOfNodes-1,a);
        for j = 1:a
            edg = perm(pos(j),:);
            set_edges = [set_edges;edg];
        end      
   else
        pos = randperm(NumberOfNodes-1,a)+((NumberOfNodes-1)*(i-1));
        for j = 1:a
            edg = perm(pos(j),:);
            set_edges = [set_edges;edg];
        end  
   end
end




%% Plotting of the Graph
 s1_30 = set_edges(:,1);
 t1_30 = set_edges(:,2);
G = digraph(s1_30,t1_30,'omitselfloops');
figure
plot(G,'Layout','force')
bins = conncomp(G);
n_bin = sum(bins);

% if n_bin = NumberofNodes then the graph is strongly connected.
 
%% Remove the unnecessary edges from the graphs. 
%  s1 = s;
% t1 = t;
% G = digraph(s1,t1,'omitselfloops');
% figure
% plot(G,'Layout','force')
% bins = conncomp(G);
% n_bin = sum(bins);

