% Erin Carson
% pcacg.m


% Run the s-step left-preconditioned CACG method to solve Ax=b using
% monomial basis

%Input:
%A: square, sparse matrix with dimension n
%b: right hand side of system to solve, Ax=b; vector of dimension n
%s: number of inner-loop iterations per outer loop; the "s" in "s-step
%methods"
%x0: initial guess for solution, vector of dimension n
%M: left preconditioner to be applied
%maxits: maximum number of iterations to complete before returning
%tol: convergence criteria for computed residual 2-norm


%Output:
%results struct stores:
%r_exact_norm: 2-norm of true residual computed in each iteration
%(results.r_exact_norm)
%r_comp_norm: 2-norm of computed residual computed in each iteration
%(results.r_comp_norm)
%x: approximate solution computed in each iteration
%(results.x)

function results = pcacg(A, b, s, x0, M, maxits, tol)

%Size of matrix
N = size(A,1);

%Set initial values for vectors
r0 = b - A*x0;
x(:,1)  = x0;
r(:,1)  = r0;
z(:,1) = M\r0;
p(:,1)  = z(:,1);

%Set outer loop iteration count to 0
k = 0;

%Set total number of iterations to 0
its = 0;

%Initialize initial true and computed residuals
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);
results.x = x0;


%Begin the iterations
while its < maxits
    
    %Break out of the loop if we have converged
    if(results.r_comp_norm(its+1) <= tol)
        break;
    end
    
    %Compute basis vectors with starting vector p
    P(:,1) = p(:,s*k+1);
    for i = 2:2*s+1
        if mod(i,2) == 0
            P(:,i) = A*P(:,i-1);
        else
            P(:,i) = M\P(:,i-1);
        end      
    end
    
    %Compute basis vectors with starting vector r
    R(:,1) = r(:,s*k+1);
    for i = 2:2*s
        if mod(i,2) == 0
            R(:,i) = M\R(:,i-1);
        else
            R(:,i) = A*R(:,i-1);
        end      
    end
    
    
    %Construct block "basis change" matrix
    Tp = diag(ones(2*s,1),-1);
    Tp = Tp(1:2*s+1,1:2*s);
    Tr = diag(ones(2*s-1,1),-1);
    Tr = Tr(1:2*s,1:2*s-1);
    Tt  = sparse([Tp, zeros(2*s+1,2*s+1); zeros(2*s,2*s+1), Tr(1:2*s,1:(2*s-1)),zeros(2*s,1)]);
    
    % Initialize CG coordinate vectors for current outer loop
    p_c = [[1;zeros(4*s,1)],zeros(4*s+1,s)];
    r_c = [[zeros(2*s+1,1);1;zeros(2*s-1,1)],zeros(4*s+1,s)];
    z_c = Tt*r_c;
    x_c = zeros(4*s+1,s+1);
    
    % Compute Gram matrix
    G = [P,R]'*[P,R];
    
    %Begin s inner iterations
    for j = 1:s
        
        if (its>=maxits)
            break;
        end
    
        %increase iteration count
        its = its + 1;
      
        %Compute scalar alpha using Gram matrix for inner products
        alpha(its) = (z_c(:,j)'*G*r_c(:,j)) / (p_c(:,j)'*G*Tt*p_c(:,j));
        
        %Update x coordinate vector
        x_c(:,j+1) = x_c(:,j) + alpha(its)*p_c(:,j);
        
        %Perform basis change to compute x vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        x(:,s*k+j+1) = [P,R]*x_c(:,j+1) + x(:,s*k+1);
        
        %Update r coordinate vector
        r_c(:,j+1) = r_c(:,j) - alpha(its)*Tt*p_c(:,j);
        
        
        %Perform basis change to compute r vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        r(:,s*k+j+1)  = [P,R]*r_c(:,j+1);
        
        %Update z coordinate vector
        z_c(:,j+1) = Tt*r_c(:,j+1);
        
        %Perform basis change to compute z vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        z(:,s*k+j+1)  = [P,R]*z_c(:,j+1);
        
        %Compute scalar beta using Gram matrix for inner products
        beta(its) = (z_c(:,j+1)'*G*r_c(:,j+1))/ (z_c(:,j)'*G*r_c(:,j));
        
        %Update p coordinate vector
        p_c(:,j+1) = z_c(:,j+1) + beta(its)*p_c(:,j);
        
        %Perform basis change to compute p vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        p(:,s*k+j+1)  = [P,R]*p_c(:,j+1);
        
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
    k=k+1;
    
end



