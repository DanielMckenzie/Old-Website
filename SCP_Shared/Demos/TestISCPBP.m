% Testing the iterated, bipartite version of single cluster pursuit. Finds 
% all clusters, as long as they are the same size!
% For further details see the function SCPBP.

clear, clc, close all
cd ../
addpath('GenerateERGraph','Functions')
warning off

% Parameters, and making the adjacency matrix
n0 = 100; %cluster row size
m0 = 200; %cluster column size
k = 10; %number of clusters
n = n0*k;
m = m0*k;
p = 0.5;
q = 20/(n-n0); % The `noise' parameters

% Creating the A matrix
[A,B] = generateBPA(n0,m0,k,p,q);
% We do not use A, the actual adjacency matrix, but note that:
% A = [0, B; B', 0]. We shall call B the incidence matrix.
I1 = mat2gray(full(B));
imshow(imcomplement(I1));
title('The incindence matrix of a G drawn from bipartite SBM. Block structure clearly visible')

% Randomly permute rows and columns of A.
index1 = 1;
pi1 = randperm(n);
pi2 = randperm(m);
B = B(pi1,pi2);
[~,pi1Inv] = sort(pi1);
[~,pi2Inv] = sort(pi2);
index1 = pi1Inv(index1);
I2 = mat2gray(full(B));
figure
imshow(imcomplement(I2));
title('The same matrix, B, but with rows and columns randomly permuted using different permutations')
C1Ltrue = pi1Inv(1:n0);
C1Rtrue = pi2Inv(1:m0);
TrueClustersL = cell(k,1);
TrueClustersR = cell(k,1);
for i = 1:k
    TrueClustersL{i} = pi1Inv((i-1)*n0 +1:i*n0);
    TrueClustersR{i} = pi2Inv((i-1)*m0 +1:i*m0);
end

 
% Run and time ISCPBP.
 tic
 [NewIndsL,NewIndsR,CL,CR]= ISCPBP(B,n0);
 ISCPBPTime = toc;
 disp(['Time taken by ISCPBP: ', num2str(ISCPBPTime)])
 B_new = B(NewIndsL, NewIndsR);
 I3 = mat2gray(full(B_new));
 figure, imshow(imcomplement(I3));
 title('Result of applying SCPBP to the scrambled B matrix. All clusters clearly visible')
 
% Compute the error
ErrorCL = 0;
ErrorCR = 0;
for i = 1:k
    tempErrorL = zeros(k,1);
    tempErrorR = zeros(k,1);
    for j = 1:k
        tempErrorL(j) = length(setdiff(TrueClustersL{i},CL{j}));
        tempErrorR(j) = length(setdiff(TrueClustersR{i},CR{j}));
    end
    realErrorL = min(tempErrorL);
    realErrorR = min(tempErrorR);
    ErrorCL = ErrorCL + realErrorL;
    ErrorCR = ErrorCR + realErrorR;
end

disp(['Number of left indices misclassified by ISCPBP ', num2str(ErrorCL)])
disp(['Number of right indices misclassified by ISCPBP ', num2str(ErrorCR)])






 