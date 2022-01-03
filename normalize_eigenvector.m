function [V1, D, v] = normalize_eigenvector(L)
%L=   [     3    -1    -1    -1; -1     2     0    -1; 0     0     1    -1; -1     0    -1     2];
[V1,D1,W1] = eig(L);
V = real(V1);
D = real(D1);
W = real(W1);
% V is the right eigenvectors. We need to normalize in such a way that the
% eigenvector corresponding to the eigenvalue 0 is 1.
for i = 1: length(V)
   for j  = 1: length(V)
    if D(i,j) <= 0.000001
        D(i,j)=0;
    end
   end
end
%here p gives the index of the column vector that corresponds to the zero
%eigenvalue. 
p =0;
for i = 1:length(D)
    if D(i,i) == 0
     p=i;
    end
end
%After we have the index value (i.e, the column of the eigenvector) we
%normalize the eigenvector.

%Following gives the sum of the values in the column vectors.
q = 0;
for i = 1:length(W)
    r1 = ((W(i,p)));
    q =q+r1;
end
%Now, dividing the matrix with the sum we get the normalized matrix
W = W/q;
W =W';
V = V/V(1,p);
u = single(V(:,p));
v = W(p,:);
V1 = V;
V1(:,p) = [];
end