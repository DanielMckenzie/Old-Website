% Testing the bipartite version of single cluster pursuit. Finds first
% cluster only.
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

 
% Run and time SCPBP.
 tic
 [C1L,C1R] = SCPBP(B,n0,index1);
 SCPBPTime = toc;
 disp(['Time taken by SCPBP: ', num2str(SCPBPTime)])
 C1Lcomp = setdiff(1:n,C1L);
 C1Rcomp = setdiff(1:m,C1R);
 newIndsL = [C1L, C1Lcomp];
 newIndsR = [C1R, C1Rcomp];
 B_new = B(newIndsL, newIndsR);
 I3 = mat2gray(full(B_new));
 figure, imshow(imcomplement(I3));
 title('Result of applying SCPBP to the scrambled B matrix. Cluster 1 is clearly visible')
 
% Computing the error
L_Error = length(setdiff(C1Ltrue,C1L));
L_Missed = length(setdiff(C1L,C1Ltrue));
disp(['Number of left indices misclassified by SCPBP ', num2str(L_Missed)])
disp(['Number of left indices missed by SCPBP ', num2str(L_Error)])

R_Error = length(setdiff(C1Rtrue,C1R));
R_Missed = length(setdiff(C1R,C1Rtrue));
disp(['Number of right indices misclassified by SCPBP ', num2str(R_Missed)])
disp(['Number of right indices missed by SCPBP ', num2str(R_Error)])

 




 