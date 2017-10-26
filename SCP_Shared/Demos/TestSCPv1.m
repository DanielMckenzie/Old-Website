% Testing Compressive Clustering using subspace pursuit.
% Test version of Single Cluster Pursuit used for numerical experiments as
% described in Lai and Mckenzie 'A Compressive Sensing Approach to
% Community Detection with Applications', 2017.
% The graph generated is drawn from the stochastic block model G(n,k,p,q).
% October 23 2017
% Be aware that SCPv3 is faster!
%

clear, clc, close all
cd ../
addpath('GenerateERGraph','Functions')
warning off

% Parameters, and making the A matrix
n = 3000;
n0 = ceil(5*sqrt(n)); %The expected size of the clusters
k = floor(n/n0); %number of clusters
p = 0.5;
q = 15/(n-n0); % The `noise' parameters

% Creating the A matrix
A = generateA(n,n0,p,q);
I1 = mat2gray(full(A)); 


% Randomly permute rows and columns of A.
index1 = 1;
pi = randperm(n);
A = A(pi,pi);
[~,piInv] = sort(pi);
index1 = piInv(index1);
Cluster1True = sort(piInv(1:n0))';
I2 = mat2gray(full(A));

% Display A before and after random permutation. Remove if timing
imshow(imcomplement(I1));
title('Adjacency matrix A of G drawn from Stochastic Block Model (SBM), with block structure clearly visible')
figure
imshow(imcomplement(I2));
title('Same matrix, but with rows and columns randomly permuted')
 
%Run and time Compressive Clustering.
 tic
 C1C = SCPv1(A,n0, index1);
 C1comp = setdiff(1:n,C1C);
 SCPv1Time = toc;
 disp(['Time taken by SCPv1 ', num2str(SCPv1Time)])
 newIndsC = [C1C, C1comp];
 A_newC = A(newIndsC, newIndsC);
 I3 = mat2gray(full(A_newC));
 figure
 imshow(imcomplement(I3));
 title('Result of appliyng ISCP to the scrambled matrix. Cluster 1 is clearly visible')
 clear newIndsC A_newC
 
% Computing the error
ErrorC = length(setdiff(Cluster1True,C1C'));
MissedC = length(setdiff(C1C',Cluster1True));
disp(['Number of indices misclassified by SCPv1 ', num2str(MissedC)])
disp(['Number of indices missed by SCPv1 ', num2str(ErrorC)])

 



 