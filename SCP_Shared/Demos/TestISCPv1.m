% Testing and timing the function ISCPv1, for finding all the clusters
% within a graph G.
% See the function `ISCPv`' for further details.

clear, clc, close all
cd ../
addpath('GenerateERGraph','Functions')
warning off

% Parameters, and making the adjacency matrix
n = 2400;
n0 = 400; %The size of the clusters
k = n/n0; %number of clusters, must be an integer!
p =0.5;
q = 35/(n-n0); % The `noise' parameter

% Creating the A matrix
A = generateA(n,n0,p,q);
I1 = mat2gray(full(A));

% Randomly permute rows and columns of A.
index1 = 1;
pi = randperm(n);
A = A(pi,pi);
I2 = mat2gray(full(A));
[~,piInv] = sort(pi);
TrueIndices = cell(k,1);
TrueClusters = cell(k,1);
for i = 1:k
    TrueClusters{i} = piInv((i-1)*n0 +1:i*n0);
    TrueIndices{i} = piInv((i-1)*n0 +1);
end

% Display A before and after random permutation. Remove if timing
imshow(imcomplement(I1));
title('Adjacency matrix A of G drawn from Stochastic Block Model (SBM), with block structure clearly visible')
figure
imshow(imcomplement(I2));
title('Same matrix, but with rows and columns randomly permuted')

% Run and time with ISCPv1
tic
[NewIndsC, Cs]=ISCPv1(A,n0);
ISCPv1Time=toc;
disp(['Time taken by ISCPv1: ', num2str(ISCPv1Time)])
A_newC = A(NewIndsC, NewIndsC);
I3 = mat2gray(full(A_newC));
clear NewIndsC A_newC

% Display unscrambled matrix
figure
imshow(imcomplement(I3));
title('Result of appliyng SCP to the scrambled matrix. All clusters are clearly visible')

%Compute the error
ErrorC = 0;
for i = 1:k
    tempError = zeros(k,1);
    for j = 1:k
        tempError(j) = length(setdiff(TrueClusters{i},Cs{j}));
    end
    realError = min(tempError);
    ErrorC = ErrorC + realError;
end
disp(['Total number of misclassified vertices: ', num2str(ErrorC)])
 