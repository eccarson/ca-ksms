% Erin Carson
% basisparams.m
% Edited 5/13/2015

%Find basis parameters to generate (s+1)-dimension Krylov basis for matrix A for selected basis type

%Input:
%A: square, sparse matrix with dimension n
%s: number of basis vectors to be generated (not counting the starting vector)
%basis_type: string denoting which basis to use. Acceptable values are
% 'monomial', 'newton', or 'chebyshev'. If something besides these strings
% entered, will default to using monomial basis. 

%Output:
%T: an (s+1)-by-s tridiagonal matrix which stores basis parameters; represents multiplication by A in
%the chosen basis
%alp, bet, gam: vectors giving the diagonal, superdiagonal, and subdiagonals of T,
%respectively, of length s, s-1, and s, respectively. 

function [alp,bet,gam,T] = basisparams( A,s,basis_type )

if (strcmp(basis_type, 'monomial'))
    alp = zeros(s,1);
    bet = zeros(s-1,1);
    gam = ones(s,1);

    %Set diagonals of T matrix
    T = zeros(s+1,s);
    T = diag(alp,0) + diag(bet,1) + diag(gam(1:end-1),-1);
    T=[T;zeros(1,s-1),gam(end)];

else

    %Compute full spectrum of A
    evec = eig(full(A));
    
    % If computing full spectrum is too computationally intensive here, the
    % above code can be commented out and the two lines below uncommented
    % instead; note however that this results in non-reproducible results from
    % run to run, since Matlab's eigs using a random starting vector.
    % Also note that this is only recommended for use when user is sure
    % that the eigenvalues of A are real
    % mx = eigs(A,1,'LM');
    % mn = eigs(A,1,'SM');
    
    %If all eigenvalues are real, call function to compute real basis
    %parameters
    if(numel(find(imag(evec)>0))==0)
        [alp,bet,gam,T] = basisparamsreal(s,evec,basis_type);
        
    %Otherwise, if there exists at least one complex eigenvalue, call
    %function to compute complex basis parameters
    else
        [alp,bet,gam,T] = basisparamscomplex (s,evec,basis_type);
    end
end

