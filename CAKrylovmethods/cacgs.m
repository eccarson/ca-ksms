% Erin Carson
% cacgs.m
% Edited 5/27/2015

% Run the s-step CACGS method to solve Ax=b for nonsymmetric A

%Input:
%A: square, sparse matrix with dimension n
%b: right hand side of system to solve, Ax=b; vector of dimension n
%s: number of inner-loop iterations per outer loop; the "s" in "s-step
%methods"
%x0: initial guess for solution, vector of dimension n
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

function results = cacgs(A, b, s, x0, maxits, tol, basis_type)

addpath('../basiscomputation/')

%Size of matrix
N = size(A,1);

%Set initial values for vectors
r0 = b - A*x0;
rt0 = r0;
x(:,1)  = x0;
r(:,1)  = r0;
p(:,1)  = r0;
u(:,1)  = r0;
delta(1) = (rt0'*r(:,1));

%Set outer loop iteration count to 0
k = 0;

%Set total number of iterations to 0
its = 0;

%Initialize initial true and computed residuals
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);

%Compute real basis parameters
[alp, bet, gam, T] = basisparams(A, 2*s, basis_type);

%Begin the iterations
while its < maxits
    
    %Break out of the loop if we have converged
    if(results.r_comp_norm(its+1) <= tol)
        break;
    end
    
    %Compute Krylov basis with starting vector p
    P = computeBasis(A,p(:,s*k+1),2*s,alp,bet,gam);
    
    %Compute Krylov basis with starting vector r
    R = computeBasis(A,r(:,s*k+1),2*s-2,alp,bet,gam);
    
    %Compute Krylov basis with starting vector q
    U = computeBasis(A,u(:,s*k+1),2*s-1,alp,bet,gam);
    
    %Construct block "basis change" matrix
    Tt  = sparse([T, zeros(2*s+1,4*s); zeros(2*s,2*s+1), T(1:end-1,1:end-1),zeros(2*s,2*s); zeros(2*s-1,4*s+1), T(1:end-2,1:end-2), zeros(2*s-1,1)]);
    
    % Initialize CGS coordinate vectors for current outer loop
    p_c = [[1;zeros(6*s-1,1)],zeros(6*s,s)];
    u_c = [[zeros(2*s+1,1);1;zeros(4*s-2,1)],zeros(6*s,s)];
    r_c = [[zeros(4*s+1,1);1;zeros(2*s-2,1)],zeros(6*s,s)];
    x_c = zeros(6*s,s+1);
    
    % Compute Gram vector
    g = rt0'*[P,U,R];
    
    %Begin s inner iterations
    for j = 1:s
        
        if (its >= maxits)
            break;
        end
    
        %increase iteration count
        its = its + 1;
      
        %Compute scalar rho using Gram vector for inner products
        alpha(its) = delta(its)/(g*Tt*p_c(:,j));
        
        q_c(:,j) = u_c(:,j) - alpha(its)*Tt*p_c(:,j);
        
        x_c(:,j+1) = x_c(:,j) + alpha(its)*(u_c(:,j)+q_c(:,j));
        x(:,s*k+j+1) = [P,U,R]*x_c(:,j+1) + x(:,s*k+1);
        
        r_c(:,j+1) = r_c(:,j) - alpha(its)*Tt*(u_c(:,j)+q_c(:,j));
        r(:,s*k+j+1)  = [P,U,R]*r_c(:,j+1);
        
        delta(its+1) = g*r_c(:,j+1);
        
        beta(its) = delta(its+1)/delta(its);
        
        u_c(:,j+1) = r_c(:,j+1) + beta(its)*q_c(:,j);
        u(:,s*k+j+1)  = [P,U,R]*u_c(:,j+1);
        
        p_c(:,j+1) = u_c(:,j+1) + beta(its)*(q_c(:,j)+ beta(its)*p_c(:,j));
        p(:,s*k+j+1)  = [P,U,R]*p_c(:,j+1);

        
        %Compute and store true residual norm (note we wouldn't do this in
        %the inner loop in practice)
        results.r_exact_norm(its+1) = norm(b-A*x(:,s*k+j+1));
        
        %Compute and store computed residual norm (note we wouldn't do this, at least in this way, in
        %the inner loop in practice)
        results.r_comp_norm(its+1) = norm(r(:,s*k+j+1));
        
        %Store current solution
        results.x = x(:,s*k+j+1);
        

    end
    
    %Increment k, outer iteration index
    k = k+1;
    
end



