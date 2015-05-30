% Erin Carson
% cadcg.m
% Edited 5/28/2015

% Run the s-step deflated CACG method to solve Ax=b

%Input:
%A is an nxn matrix, b is RHS (Ax=b)
%b is RHS
%s is the basis size (number of times we've unrolled the loop)
%t is the max number of outer iterations (s*t = max iterations)
%x_g is the initial guess (zero vector is ok here)
%kk is number of deflation vectors
%tol is convergence requirement (r norm)
%basis is string that is either 'monomial', 'newton', or 'chebyshev'
%hasW is 1 if deflation vectors provided
%Win is deflation vectors

%Output:
%results struct stores:
%r_exact_norm: 2-norm of true residual computed in each iteration
%(results.r_exact_norm)
%r_comp_norm: 2-norm of computed residual computed in each iteration
%(results.r_comp_norm)
%x: approximate solution computed in each iteration
%(results.x)

function results = cadcg(A, b, s, maxits, tol, basis, kk, hasW, Win)

addpath('../basiscomputation/')

%If no deflation vectors are provided, use Matlab eigs to compute
%eigenvector estimates for kk smallest eigenvalues
if(~hasW)
    [W,~]= eigs(A, kk, 'SM');
%If deflation vectors are given, use the provided W
else
    W=Win;
end

%Precompute and store W'AW (could also prefactorize in practice)
waw = (W'*A*W);

%Deflate initial approximation
x0 = W*((waw)\W'*b);

%Initialize results struct
results.r_exact_norm = 0;

%Set initial values for vectors
r0 = b - A*x0;
mu0 = waw\(W'*A*r0);
p0 = r0-W*mu0;
r_2norm(1,1) = norm(r0,2);
results.x = x0;
x(:,1)  = x0;
r(:,1)  = r0;
p(:,1)  = p0;
k = 0;

%Initialize initial true and computed residuals
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);

%Compute basis parameters 
[alps2,bets2,gams2,Tp] = basisparams(A,s+1,basis);
Tr = Tp(1:end-1,1:end-1);
Tp = [Tp,zeros(s+2,1)]; Tr = [Tr,zeros(s+1,1)]; 
alps1 = alps2(1:end-1); alps = alps1(1:end-1);
bets1 = bets2(1:end-1); bets = bets1(1:end-1);
gams1 = gams2(1:end-1); gams = gams1(1:end-1);

%Precompute Krylov basis with deflation vectors as RHSs 
if(kk>0)
    [WB] = computeBasis(A,W,s-1,alps,bets,gams);
    TW = sparse(kron(Tr(1:end-1,1:end-1),eye(kk)));
else
    WB=[];
    TW=[];
end
its = 0;

%Begin the iterations
while its < maxits

    %Break out of the loop if we have converged
    if(results.r_comp_norm(its+1) <= tol)
        break;
    end
    
    %Compute Krylov basis with starting vector p
    P = computeBasis(A,p(:,s*k+1),s+1,alps2, bets2, gams2);
    
    %Compute Krylov basis with starting vector r
    R = computeBasis(A,r(:,s*k+1),s,alps1, bets1, gams1);

    %Construct T matrix
    T1  = sparse ( [Tp, zeros(s+2,s+1); zeros(s+1,s+2), Tr] );
    T = sparse([T1, zeros(2*s+3,(s)*kk); zeros((s)*kk,2*s+3), TW]);
 
    % Initialize coefficient vectors
    p_c = [[1;zeros(2*s+2+(s)*kk,1)],zeros(2*s+3+(s)*kk,s)];   
    r_c = [[zeros(s+2,1);1;zeros(s,1);zeros((s)*kk,1)],zeros(2*s+3+(s)*kk,s)];
    x_c = zeros(2*s+3+(s)*kk,s+1);
    
    %Compute Gram matrix
    V = [P,R,WB];
    G = V'*V;

    %Compute Z (Could extract from already computed G using Z = G(2*s+3+1:2*s+3+kk,:);
    Z=W'*V;
    
    %Inner iterations
    for j = 1:s
        
        if (its >= maxits)
            break;
        end
       
        %Increment iteration count
        its = its + 1;
  
        %Compute scalar alpha using Gram matrix
        alpha(its) = (r_c(:,j)'*G*r_c(:,j)) / (p_c(:,j)'*G*T*p_c(:,j)); 
        
        %Update x coordinate vector
        x_c(:,j+1) = x_c(:,j) + alpha(its)*p_c(:,j);
        
        %Perform basis change to compute x vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        x(:,s*k+j+1)  = x(:,s*k+1) + V*x_c(:,j+1);
        
        %Update r coordinate vector
        r_c(:,j+1) = r_c(:,j) - alpha(its)*T*p_c(:,j);
        
        %Perform basis change to compute r vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        r(:,s*k+j+1)  = V*r_c(:,j+1);
        
        %Compute scalar beta using Gram matrix
        beta(its) = (r_c(:,j+1)'*G*r_c(:,j+1)) / (r_c(:,j)'*G*r_c(:,j));

        %Compute mu using precomputed W'AW
        muj = (waw)\(Z*T*r_c(:,j+1));
        
        %Update p coordinate vector
        p_c(:,j+1) = r_c(:,j+1) + beta(its)*p_c(:,j) - [zeros(2*s+3,1);muj;zeros((s-1)*kk,1)];
      
        %Perform basis change to compute r vector in standard basis (note we wouldn't need to do this
        %in the inner loop in practice)
        p(:,s*k+j+1)  = V*p_c(:,j+1);
        
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


