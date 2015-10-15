% Erin Carson
% calsqr_alt.m
% Edited 7/8/2015

% Run the s-step CALSQR method to solve min||Ax=b|| 

%Input:
%A: mxn sparse matrix 
%b: right hand side of least squares system to solve, Ax=b; vector of dimension n
%s: number of inner-loop iterations per outer loop; the "s" in "s-step
%methods"
%maxits: maximum number of iterations to complete before returning
%tol: convergence criteria for computed residual 2-norm
%basis_type: string denoting which basis to use. Acceptable values are
% 'monomial', 'newton', or 'chebyshev'. If something besides these strings
% entered, will default to using monomial basis. 

%Output:
%results struct stores:
%r_exact_norm: 2-norm of true residual computed in each iteration
%(results.r_exact_norm)
%r_comp_norm: 2-norm of computed residual computed in each iteration
%(results.r_comp_norm)
%x: approximate solution computed in each iteration
%(results.x)

function [results] = calsqr_alt(A, b, s, maxits, tol, basis_type)

addpath('../basiscomputation/')

%Size of matrix
mm = size(A,1);
n = size(A,2);

%Initialize data structures
U=zeros(mm,maxits);
V=zeros(n,maxits);
W=zeros(n,maxits);
X=zeros(n,maxits);

%Initialize values
utemp = b;
beta(1) = norm(b);
U(:,1) = utemp/beta(1);
vtemp = A.'*U(:,1);
alpha(1) = norm(vtemp);
V(:,1)= vtemp/alpha(1);

%Initialize LSQR scalars
phit(1) = beta(1);
rhot(1) = alpha(1);
W(:,1) = V(:,1);
X(:,1) = zeros(n,1);

%Initialize initial true residuals
results.r_exact_norm(1) = norm(b-A*X(:,1));

%Set outer loop iteration count to 0
k = 0;

%Set total number of iterations to 0
m = 0;

%Compute real basis parameters
[alp, bet, gam, T] = basisparams(A, 2*s, basis_type);
%evec=eig(full(A));
%[alp,bet,gam,T] = basisparamscomplex (2*s,evec,basis_type);
T = [T, zeros(2*s+1,1)];

%Begin the iterations!
while m < maxits
    
    %Break out of the loop if we have converged
    if(results.r_exact_norm(m+1) <= tol)
        break;
    end
    
    %Compute alternating bases with starting vectors v and u
    [Z,Y] = computeBasis_alt(A,V(:,m+1),U(:,m+1),2*s,alp,bet,gam);
    %Y = [Y(:,1:end-1),zeros(n,1)];
    
    %Initialize coordinate vectors for current outer loop
    Itmp = eye(2*s+1); 
    vcoeff = Itmp(:,1);
    ucoeff = Itmp(:,1);
    wcoeff = [zeros(2*s+1,1);1];
    xcoeff = zeros(2*s+2,1);   

    %Compute Gram matrix
    Gy = Y'*Y;
    Gz = Z'*Z;
    
    %Store current x and w vectors
    XM = X(:,m+1);
    WM = W(:,m+1);
    
    %Begin s inner iterations
    for j = 1:s
           
        if (m >= maxits)
            break;
        end
    
        %increase iteration count
        m = m + 1;
 
        %Update unnormalized u coordinate vector
        utmp = T*vcoeff(:,j)-alpha(m)*ucoeff(:,j);
    
        %Compute scalar beta using Gram matrix for inner products
        beta(m+1) = sqrt(utmp'*Gy*utmp);
        
        %Normalize u coordinate vector
        ucoeff(:,j+1) = utmp/beta(m+1);
       
        %Perform basis change to compute u vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        U(:,m+1) = Y*ucoeff(:,j+1);
        
        %Update unnormalized v coordinate vector
        vtmp = T*ucoeff(:,j+1) - beta(m+1)*vcoeff(:,j);
        
        %Compute scalar alpha using Gram matrix for inner products
        alpha(m+1) = sqrt(vtmp'*Gz*vtmp);
        
        %Normalize v coordinate vector
        vcoeff(:,j+1) = vtmp/alpha(m+1);
        
        %Perform basis change to compute v vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        V(:,m+1) = Z*vcoeff(:,j+1);
        
        %Update LSQR scalars
        rho(m) = (rhot(m)^(2) + beta(m+1)^2)^(1/2);
        c(m) = rhot(m)/rho(m);
        ss(m) = beta(m+1)/rho(m);
        theta(m+1) = ss(m)*alpha(m+1);
        rhot(m+1) = -c(m)*alpha(m+1);
        phi(m) = c(m)*phit(m);
        phit(m+1) = ss(m)*phit(m);
        
        %Update solution coordinate vector
        xcoeff(:,j+1) = xcoeff(:,j) + phi(m)/rho(m)*wcoeff(:,j);
        
        %Perform basis change to compute x vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        X(:,m+1) = [Z,WM]*xcoeff(:,j+1) + XM;
        
        %Update w coordinate vector
        wcoeff(:,j+1) = [vcoeff(:,j+1);0] - (theta(m+1)/rho(m))*wcoeff(:,j);
        
        %Perform basis change to compute w vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        W(:,m+1) = [Z,WM]*wcoeff(:,j+1);
       
        %Compute and store true residual norm (note we wouldn't do this in
        %the inner loop in practice)
        results.r_exact_norm(m+1) = norm(b-A*X(:,m+1));
        
        %Compute and store computed residual norm (note we wouldn't do this, at least in this way, in
        %the inner loop in practice)
        results.r_comp_norm(m) = phit(m+1);
        
        %Store current solution
        results.x = X(:,m+1);
        
    end
 
end

end
