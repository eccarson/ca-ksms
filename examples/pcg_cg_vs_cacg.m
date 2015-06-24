% Erin Carson
% pcg_cg_vs_cacg.m
% Edited 6/23/2015

addpath('../CAKrylovmethods/')
addpath('../utils/')

% Create poisson problem matrix (2D 5-pt stencil)
n = 16;
A = gallery('poisson',n);
N=n^2;

% Set the right hand side b
x_true = ones(N,1)./sqrt(N);
b = A*x_true;

% Set s, maximum number of iterations, and tolerance
maxits = N;
pcmaxits = 10;
tol = 1e-14;
pctol = 1e-16;
s=20;
basis='monomial';

% Set initial solution to 0 vector
x0 = zeros(N,1);

% Call PCG_CG method
results_pcg_cg = pcg_cg(A, b, x0, maxits, tol, pcmaxits, pctol);

% Call Matlab's CG method for comparison
results_cg = cg(A, b, x0, maxits, tol);

% Call PCG_CACG method
results_pcg_cacg = pcg_cacg(A, b, s, x0, maxits, tol, pcmaxits, pctol, basis);

% Generate plot showing convergence of computed residual
figure();
semilogy(1:numel(results_pcg_cg.r_comp_norm), results_pcg_cg.r_comp_norm,'b:');
hold on;
semilogy(1:numel(results_pcg_cacg.r_comp_norm), results_pcg_cacg.r_comp_norm,'r-');
hold on;
semilogy(1:numel(results_cg.r_comp_norm), results_cg.r_comp_norm,'k.');
title('CG Convergence');
xlabel('Iteration');
ylabel('Residual 2-norm');


