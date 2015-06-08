% Erin Carson
% cgs.m
% Edited 1/14/2015

% Run the CGS method to solve Ax=b

%Input:
%A: square, sparse matrix with dimension n
%b: right hand side of system to solve, Ax=b; vector of dimension n
%x0: initial guess for solution, vector of dimension n
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

function results = cgs(A, b, x0, maxits, tol)

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

%Set total number of iterations to 0
its = 0;

%Initialize initial true and computed residuals
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);


%Begin the iterations
while its < maxits
    
    %Break out of the loop if we have converged
    if(results.r_comp_norm(its+1) <= tol)
        break;
    end
     
    %increase iteration count
    its = its + 1;

    alpha(its) = delta(its)/(rt0'*A*p(:,its));

    q(:,its) = u(:,its) - alpha(its)*A*p(:,its);

    x(:,its+1) = x(:,its) + alpha(its)*(u(:,its)+q(:,its));

    r(:,its+1) = r(:,its) - alpha(its)*A*(u(:,its)+q(:,its));

    delta(its+1) = rt0'*r(:,its+1);

    beta(its) = delta(its+1)/delta(its);

    u(:,its+1) = r(:,its+1) + beta(its)*q(:,its);

    p(:,its+1) = u(:,its+1) + beta(its)*(q(:,its)+ beta(its)*p(:,its));


    %Compute and store true residual norm 
    results.r_exact_norm(its+1) = norm(b-A*x(:,its+1));

    %Compute and store computed residual norm 
    results.r_comp_norm(its+1) = norm(r(:,its+1));

    %Store current solution
    results.x = x(:,its+1);
          
end



