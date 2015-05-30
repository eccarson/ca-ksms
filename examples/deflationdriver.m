% Erin Carson
% deflationdriver.m
% Edited 5/30/2015

%A is system matrix (SPD)
%b is right hand side to use for all systems
%basis is either 'monomial', 'newton', or 'chebyshev'
%s is #iterations per outer loop
%xlim is maximum #iterations

%example call: deflationdriver(gallery('poisson',16), rand(256,1), 'monomial', 8, 100)

function deflationdriver(A, b, basis, s, xlim)

addpath('../CAKrylovmethods/')
addpath('../utils/')

tol=1e-9;

%max deflation vectors to use
nn=10;

%get spectrum of A
[V,D] = eig(full(A));
[~,I] = sort(abs(diag(D)));

%get nn smallest eigenvectors of A
V = V(:,I(1:nn));

%Set W's for 1, 2, and 10 deflation vectors
W1 = V(:,1:1);
W2 = V(:,1:2);
W10 = V(:,1:10);

%No deflation, CG
results = cadcg(A, b, 1, xlim, tol, basis, 0, 0, 0);
resultsN = results.r_comp_norm;

%Deflation, CG
results = cadcg(A, b, 1, xlim, tol, basis, size(W1,2), 1, W1);
resultsN1 = results.r_comp_norm;

results = cadcg(A, b, 1, xlim, tol, basis, size(W2,2), 1, W2);
resultsN2 = results.r_comp_norm;

results = cadcg(A, b, 1, xlim, tol, basis, size(W10,2), 1, W10);
resultsN3 = results.r_comp_norm;

%No deflation, CACG
results = cadcg(A, b, s, xlim, tol, basis, 0, 0, 0);
resultsNC = results.r_comp_norm;

%Deflation, CACG
results = cadcg(A, b, s, xlim, tol, basis, size(W1,2), 1, W1);
resultsC1 = results.r_comp_norm;

results = cadcg(A, b, s, xlim, tol, basis, size(W2,2), 1, W2);
resultsC2 = results.r_comp_norm;

results = cadcg(A, b, s, xlim,  tol, basis, size(W10,2),1,W10);
resultsC3 = results.r_comp_norm;


%Plot results
sz=max([numel(resultsN), numel(resultsNC), numel(resultsC1), numel(resultsC2), numel(resultsC3),...
    numel(resultsN1), numel(resultsN2), numel(resultsN3)]);
BB=zeros(sz,10);
BB(:,1) = 1:sz;
BB(1:numel(resultsN),2)=resultsN./resultsN(1);
BB(1:numel(resultsNC),3)=resultsNC./resultsNC(1);
BB(1:numel(resultsC1),4)=resultsC1./resultsC1(1);
BB(1:numel(resultsC2),5)=resultsC2./resultsC2(1);
BB(1:numel(resultsC3),6)=resultsC3./resultsC3(1);
BB(1:numel(resultsN1),7)=resultsN1./resultsN1(1);
BB(1:numel(resultsN2),8)=resultsN2./resultsN2(1);
BB(1:numel(resultsN3),9)=resultsN3./resultsN3(1);


for i = 2:9
    ind = find(BB(:,i)<tol,1);
    if(numel(ind)==1)
        BB(ind+1:end,i)=0;
    else
       BB(end-1:end,i)=0;
    end
end


figure()
axes1 = axes('YScale','log','YMinorTick','on',...
    'FontSize',20,...
    'FontName','Helvetica');

box(axes1,'on');
hold(axes1,'all');

s1 = semilogy(BB(:,1),BB(:,2:9),'LineWidth',1);
set(s1(1),'LineStyle','--',...
   'Color',[0 0 0],...
   'DisplayName','CG',...
   'MarkerSize',6);
set(s1(2),'LineStyle','-',...
   'Color',[0 0 0],...
   'DisplayName','CA-CG',...
   'MarkerSize',6);
set(s1(3),'LineStyle','-',...
   'Color',[1 0 0],...
   'DisplayName','CA-CG+D (c=1)',...
   'MarkerSize',6);
set(s1(4),'LineStyle','-',...
   'Color',[0 1 0],...
   'DisplayName','CA-CG+D (c=2)',...
   'MarkerSize',6);
set(s1(5),'LineStyle','-',...
   'Color',[0 0 1],...
   'DisplayName','CA-CG+D (c=3)',...
   'MarkerSize',6);
set(s1(6),'LineStyle','--',...
   'Color',[1 0 0],...
   'DisplayName','CG+D (c=1)',...
   'MarkerSize',6);
set(s1(7),'LineStyle','--',...
   'Color',[0 1 0],...
   'DisplayName','CG+D (c=2)',...
   'MarkerSize',6);
set(s1(8),'LineStyle','--',...
   'Color',[0 0 1],...
   'DisplayName','CG+D (c=3)',...
   'MarkerSize',6);

% Create xlabel
xlabel('Iteration','FontSize',20,'FontName','Helvetica');

% Create ylabel
ylabel('Residual (2-norm)','FontSize',20,'FontName','Helvetica');

axis([0 xlim 1e-10 1e1])

end