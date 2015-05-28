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
p(:,1)  = zeros(N,1);
q(:,1)  = zeros(N,1);
rho_old = 1;

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
    R = computeBasis(A,r(:,s*k+1),2*s,alp,bet,gam);
    
    %Compute Krylov basis with starting vector q
    Q = computeBasis(A',q(:,s*k+1),2*s,alp,bet,gam);
    
    %Construct block "basis change" matrix
    Tt  = sparse([T, zeros(2*s+1,4*s+3); zeros(2*s+1,2*s+1), T,zeros(2*s+1,2*s+2); zeros(2*s+1,4*s+2), T, zeros(2*s+1,1)]);
    
    % Initialize CGS coordinate vectors for current outer loop
    p_c = [[1;zeros(6*s+2,1)],zeros(6*s+3,s)];
    r_c = [[zeros(2*s+1,1);1;zeros(4*s+1,1)],zeros(6*s+3,s)];
    q_c = [[zeros(4*s+2,1);1;zeros(2*s,1)],zeros(6*s+3,s)];
    x_c = zeros(6*s+3,s+1);
    
    % Compute Gram vector
    g = rt0'*[P,R,Q];
    
    %Begin s inner iterations
    for j = 1:s
        
        if (its >= maxits)
            break;
        end
    
        %increase iteration count
        its = its + 1;
      
        %Compute scalar rho using Gram vector for inner products
        rho_new = g*r_c(:,j);
        
        %Compute scalar beta
        beta(its) = rho_new/rho_old;
        
        %Update u (auxiliary) coordinate vector
        u_c(:,j) = r_c(:,j) + beta(its)*q_c(:,j);
        
        %Update p coordinate vector
        p_c(:,j+1) = u_c(:,j) + beta(its)*(q_c(:,j)+beta(its)*p_c(:,j));
        
        %Perform basis change to compute p vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        p(:,s*k+j+1)  = [P,R,Q]*p_c(:,j+1);
                
        %Update v (auxiliary) coordinate vector       
        v_c(:,j) = Tt*p_c(:,j+1);

        %Compute scalar sigma using Gram vector for inner products
        sigma(its) = g*v_c(:,j);
        
        %Compute scalar alpha
        alpha(its) = rho_new / sigma(its);
        
        %Update q coordinate vector
        q_c(:,j+1) = u_c(:,j) - alpha(its)*v_c(:,j);
        
        %Perform basis change to compute q vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        q(:,s*k+j+1)  = [P,R,Q]*q_c(:,j+1);
        
        %Update r coordinate vector
        r_c(:,j+1) = r_c(:,j) - alpha(its)*Tt*(u_c(:,j) + q_c(:,j+1));
        
        %Perform basis change to compute r vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        r(:,s*k+j+1)  = [P,R,Q]*r_c(:,j+1);

        %Update x coordinate vector
        x_c(:,j+1) = x_c(:,j) + alpha(its) *(u_c(:,j) + q_c(:,j+1));

        %Perform basis change to compute x vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        x(:,s*k+j+1) = [P,R,Q]*x_c(:,j+1) + x(:,s*k+1);
        
        %Compute and store true residual norm (note we wouldn't do this in
        %the inner loop in practice)
        results.r_exact_norm(its+1) = norm(b-A*x(:,s*k+j+1));
        
        %Compute and store computed residual norm (note we wouldn't do this, at least in this way, in
        %the inner loop in practice)
        results.r_comp_norm(its+1) = norm(r(:,s*k+j+1));
        
        %Store current solution
        results.x = x(:,s*k+j+1);
        
        %Store rho for next iteration
        rho_old = rho_new;
        
    end
    
    %Increment k, outer iteration index
    k = k+1;
    
end



