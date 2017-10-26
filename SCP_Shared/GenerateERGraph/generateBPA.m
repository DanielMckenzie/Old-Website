function [A,B] = generateBPA(n0,m0,k,p,q)
% This script generates the adjacency matrix of a bipartite graph G drawn
% from the Bipartite Stochastic Block Model G(n0k,m0k,k,p,q). 
% k is the number of clusters and each cluster is of size n0 by m0.

B = zeros(n0*k,m0*k);

for j = 1:k
    Btemp = +(rand(n0,m0) <= p);
    Bnoise1 = +(rand((k-j)*n0,m0) <=q);
    Bnoise2 = +(rand(n0,(k-j)*m0) <=q);
    B((j-1)*n0+1:j*n0,(j-1)*m0+1:j*m0) = Btemp;
    B(j*n0+1:end,(j-1)*m0+1:j*m0) = Bnoise1;
    B((j-1)*n0+1:j*n0,j*m0+1:end) = Bnoise2;
end
A = [zeros(k*n0,k*n0), B; B', zeros(k*m0,k*m0)];
end