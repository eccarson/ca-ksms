% Erin Carson
% basisparams.m
% Edited 1/14/2015

% Compute s-step basis parameters for basis specified in basis_type, based on
% properties of A

%Input:
%A: matrix we want to compute basis parameters for
%s: number of basis parameters to compute
%basis_type: string denoting which basis to use. Acceptable values are
% 'monomial', 'newton', or 'chebyshev'. If something besides these strings
% entered, will default to using monomial basis. 

%Output:
%alp, bet, gam: vectors which form diagonal, superdiagonal, and subdiagonal
%of T, respectively
%T: matrix of size (s+1) by s; represents multiplication by A in chosen
%Krylov basis

function [alp,bet,gam,T] = basisparamsreal(s,evec,basis_type)

%Get largest and smallest eigenvalues of A
mx = max(real(evec));
mn = min(real(evec));

%Compute Newton basis
if(strcmp(basis_type, 'newton'))
    
    %Find leja points depending on extremal ritz values
    lejax = lejapointsreal(s, mn, mx);
    
    %Set basis parameters (diagonals of tridiagonal basis change matrix)
    alp = lejax;
    bet = zeros(s-1,1);
    gam = ones(s,1);
    
%Compute Chebyshev basis
elseif(strcmp(basis_type, 'chebyshev'))
    
    %Compute ellipse parameters
    aa = (mx-mn);
    cc = (mx+mn)/2;
    bb = 0;
    dd = sqrt((mx-cc)^2-bb^2);
    gg = max((mx-cc),bb);
    
    %Set basis parameters (diagonals of tridiagonal basis change matrix)
    alp = cc.*ones(s,1);
    bet = (dd)^2/(4*gg).*ones(s-1,1);
    gam = [2*gg; gg.*ones(s-1,1)];
    
%Default to monomial basis
else
    
    %Set basis parameters (diagonals of tridiagonal basis change matrix)
    alp = zeros(s,1);
    bet = zeros(s-1,1);
    gam = ones(s,1);
    
end

%Set diagonals of T matrix
T = zeros(s+1,s);
T = diag(alp,0) + diag(bet,1) + diag(gam(1:end-1),-1);
T=[T;zeros(1,s-1),gam(end)];


