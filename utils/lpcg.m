% Erin Carson
% lpcg.m
% Edited 2021

% Run the left preconditioned CG method to solve Ax=b

%Input:
%A: square, sparse matrix with dimension n
%b: right hand side of system to solve, Ax=b; vector of dimension n
%x0: initial guess for solution, vector of dimension n
%M: left preconditioner
%maxits: maximum number of iterations to complete before returning; should
%be a multiple of s
%tol: convergence criteria for computed residual 2-norm

%Output:
%results struct stores:
%r_exact_norm: 2-norm of true residual computed in each iteration
%(results.r_exact_norm)
%r_comp_norm: 2-norm of computed residual computed in each iteration
%(results.r_comp_norm)
%x: approximate solution computed in each iteration
%(results.x)

function results = cg(A, b, x0, M, maxits, tol)

%Size of matrix
N = size(A,1);

%Set initial values for vectors
r0 = b - A*x0;
x(:,1)  = x0;
r(:,1)  = r0;
z(:,1) = M\r0;
p(:,1)  = z(:,1);

%Set total number of iterations to 0
its = 0;

%Initialize initial true and computed residuals
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);
results.x=x0;


%Begin the iterations
while its < maxits
    
    %Break out of the loop if we have converged
    if(results.r_comp_norm(its+1) <= tol)
        break;
    end
     
    %increase iteration count
    its = its + 1;
       
    %Compute scalar alpha
    alpha(its) = z(:,its)'*r(:,its)/(p(:,its)'*A*p(:,its));

    %Update x vector
    x(:,its+1) = x(:,its) + alpha(its)*p(:,its);

    %Update r vector
    r(:,its+1) = r(:,its) - alpha(its)*A*p(:,its);

    %Update preconditioned residual
    z(:,its+1) = M\r(:,its+1);
    
    %Compute scalar beta
    beta(its) = (z(:,its+1)'*r(:,its+1))/ (z(:,its)'*r(:,its));

    %Update p coordinate vector
    p(:,its+1) = z(:,its+1) + beta(its)*p(:,its);

    %Compute and store true residual norm 
    results.r_exact_norm(its+1) = norm(b-A*x(:,its+1));

    %Compute and store computed residual norm 
    results.r_comp_norm(its+1) = norm(r(:,its+1));

    %Store current solution
    results.x = x(:,its+1);
          
end



