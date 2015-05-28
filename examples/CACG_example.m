% Erin Carson
% CACG_example.m
% Edited 5/13/2015

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
s = 4;
maxits = N;
tol = 1e-16;

% Set initial solution to 0 vector
x0 = zeros(N,1);

% Set basis_info.type to either 'monomial', 'newton', or 'chebyshev',
% depending on what basis you want to use.
basis_type = 'monomial';

% Call CACG method
results = cacg(A, b, s, x0, maxits, tol, basis_type);

% Call Matlab's CG method for comparison
results_cg = cg(A, b, x0, maxits, tol);

% Generate plot showing convergence of computed residual
figure();
semilogy(1:numel(results.r_comp_norm), results.r_comp_norm,'-r');
hold on;
semilogy(1:numel(results_cg.r_comp_norm), results_cg.r_comp_norm,'-k');
title(strcat('CACG Convergence, s = ',num2str(s),' , basis = ',basis_type));
xlabel('Iteration');
ylabel('Residual 2-norm');


