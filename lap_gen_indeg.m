% This function generates the indegree laplacian of the graph.
% The inputs are :
% N = number of nodes
% s = head of the edge
% t = tail of the edge
% The output L is the Laplacian matrix.
function [L] = lap_gen_indeg (N,s,t)
% example of the input.
% N=8;
% s = [ 1     1     1     2     2     3     3     4     5     5     6     6     6     6     7     7     7 7     8     8];
% t = [3     4     5     6     1     7     2     8     3     1     3     1     5     4     2     8     4   5     6     3];
%% Adjacency Matrix
% This generates the adjancecy matirx. 
adj = zeros(N,N);
for i = 1:length(s)
   adj(t(i),s(i)) =1;
end
%% Degree Matrix
% This generates the degree matrix.
deg = zeros(N,N);
dg = 0;
for i = 1:N
    for j  = 1:N
       dg = dg+ adj(i,j);
    end
    deg(i,i) = dg;
    dg = 0;
end

%% The Laplacian is given by the following
L = deg-adj;
end
