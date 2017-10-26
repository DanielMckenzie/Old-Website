function A = generateAUnequalSize(n,nvec,p,q)
% This script generates the adjacency matrix of a graph G drawn from the 
% Stochastic Block Model G(n,nvec,p,q). where nvec is a length k vector of
% the sizes of the clusters.

k = length(nvec);
A = sparse(n,n);
A(1:nvec(1),1:nvec(1)) = randsym(nvec(1),p);
B = +(sprand(nvec(1),n-nvec(1),q)>0);
A(1:nvec(1), nvec(1)+1:end) = B;
A(nvec(1)+1:end,1:nvec(1)) = B'; 

for i = 2:k
    Start = sum(nvec(1:i-1))+1;
    Finish = sum(nvec(1:i));
    A(Start:Finish,Start:Finish) = randsym(nvec(i),p);
    B = +(sprand(nvec(i),n - Finish,q)>0); % The 'noise matrix'
    A(Start:Finish,Finish+1:end) = B;
    A(Finish+1:end,Start:Finish) = B';  
end
end