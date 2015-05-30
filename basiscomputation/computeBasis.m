% Erin Carson
% CACG_example.m
% Edited 1/14/2015

%Input:
%A: matrix with which to compute Krylov basis
%x: starting vector
%s: number of SpMVs to perform; resulting basis is of size s+1
%alp, bet, gam: parameters for basis recurrence

%Output:
%V: matrix of dimension n by (s+1) whose columns are basis for the desired
%Krylov subspace

function V = computeBasis(A,x,s,alp,bet,gam)

%dimension of matrix A
n = size(A,1);

%columns in RHS
r = size(x,2);

%Set initial vector
V(:,1:r) = x;

if(s>0)
    %Compute first basis vector    
    V(:,r+1:2*r) = (1/gam(1))*(A-alp(1)*eye(n))*V(:,1:r);

    %Compute remaining s-1 basis vectors
    for j = 2:s
        V(:,((j)*r+1):((j+1)*r)) = (1/gam(j))*((A-alp(j)*eye(n))*V(:,((j-1)*r+1):((j)*r)) - bet(j-1)*V(:,((j-2)*r+1):((j-1)*r)));
    end
end


end