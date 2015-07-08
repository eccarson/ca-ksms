% Erin Carson
% computeBasis_alt.m
% Edited 7/8/2015

%Input:
%A: matrix with which to compute Krylov basis
%z1: starting vector
%y1: starting vector
%s: number of total SpMVs to perform; resulting basis is of size s+1
%alp, bet, gam: parameters for basis recurrence

%Output:
%Z,Y: matrices of dimension n by (s+1) whose columns are basis for the desired
%Krylov subspace

function [Z,Y] =  computeBasis_alt(A,z1,y1,s,alp,bet,gam)

%dimension of matrix A
n = size(A,1);

%Set initial vectors
Y(:,1) = y1;
Z(:,1) = z1;

if(s>0)
    
    %Compute first basis vector
    Z(:,2) = (A'*Y(:,1) -alp(1).*Z(:,1))./gam(1);
    Y(:,2) = (A*Z(:,1) -alp(1).*Y(:,1))./gam(1);

    %Compute remaining basis vectors
    for qq = 2:(s)
            Z(:,qq+1) = (A'*Y(:,qq) -alp(qq).*Z(:,qq) - bet(qq-1)*Z(:,qq-1))./gam(qq);
            Y(:,qq+1) = (A*Z(:,qq) -alp(qq).*Y(:,qq) - bet(qq-1)*Y(:,qq-1))./gam(qq);
    end
    
end
   
    
end