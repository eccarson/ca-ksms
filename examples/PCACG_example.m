% Erin Carson
% PCACG_example.m
% Edited 2021

addpath('../CAKrylovmethods/')
addpath('../utils/')

% Create poisson problem matrix (2D 5-pt stencil)
n = 16;
A = gallery('poisson',n);

% Set dimension of A, N = n^2
N = size(A,1);

% Set the right hand side b
x_true = ones(N,1)./sqrt(N);
b = A*x_true;

% Set s, maximum number of iterations, and tolerance
s = 8;
maxits = N;
tol = 1e-16;

% Set initial solution to 0 vector
x0 = zeros(N,1);

% Set basis_info.type to either 'monomial', 'newton', or 'chebyshev',
% depending on what basis you want to use.
basis_type = 'monomial';

% Call CACG method
results = cacg(A, b, s, x0, maxits, tol, basis_type);

% Call PCACG method with incomplete LU preconditioner
[L,U] = ilu(A);
M = L*U;
resultsp = pcacg(A, b, s, x0, M, maxits, tol);

% Call classical CG method for comparison
results_cg = cg(A, b, x0, maxits, tol);

% Call preconditioned CG method for comparison
results_cgp = lpcg(A, b, x0, M, maxits, tol);

% Generate plot showing convergence of exact residual
figure();
semilogy(1:numel(results.r_exact_norm), results.r_exact_norm,'-r');
hold on;
semilogy(1:numel(resultsp.r_exact_norm), resultsp.r_exact_norm,'--r');
semilogy(1:numel(results_cg.r_exact_norm), results_cg.r_exact_norm,'-k');
semilogy(1:numel(results_cgp.r_exact_norm), results_cgp.r_exact_norm,'--k');
title(strcat('CACG Convergence, s = ',num2str(s),' , basis = ',basis_type));
xlabel('Iteration');
ylabel('Residual 2-norm');
legend('CACG','LP-CACG', 'CG', 'LP-CG')


