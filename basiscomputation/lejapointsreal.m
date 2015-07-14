% Erin Carson
% lejapoints.m
% Edited 10/28/2012

% Compute n leja points in the real interval [a,b]

function xleja = lejapointsreal(n,a,b)

options = [];
xleja = zeros(n,1);
xmaxloc = zeros(n-2,1);
fmaxloc = xmaxloc;

xleja(1) = b;
xleja(2) = a;

xsort = xleja;
for k = 3:1:n
    xsort = sort(xsort(1:k-1));
    for i = 1:(k-2)/1
        [xmaxloc(i),fmaxloc(i)] = fminbnd(@fleja,xsort(i),xsort(i+1),options,xleja,k-1);
    end
    [fmax index] = max(-fmaxloc(1:k-2));
    xleja(k) = xmaxloc(index);
    xleja(k+1) = -xmaxloc(index);
    xsort(k) = xleja(k);
    xsort(k+1) = xleja(k+1);
end
xleja = xleja(1:n);


function z = fleja(x,xleja,n)
% the opposite of the Leja function to maximize
z = -abs(prod(x-xleja(1:n)));


function z = lejaorder(x,n)
% Leja order of the first n components of the vector x
m = length(x);
n = min(n,m);
z = zeros(size(x));
[v index] = max(abs(x));
z(1) = x(index);
temp = abs(x-z(1));
[v index] = max(temp);
z(2) = x(index);
for k = 2:n-1
    for i = 1:m
        temp(i) = temp(i)*abs(x(i)-z(k));
    end
    [v index] = max(temp);
    z(k+1) = x(index);
end
z = z(1:n);

function [min minval] = fminbnd(func,lb,ub, options, varargin)
delta = 1e-17;
gr = (sqrt(5)-1)/2;
width = (ub-lb);
out = [ lb:(width/3):ub ];
out(2) = out(4)-gr*width;
out(3) = out(1)+gr*width;
upper = feval(func,out(3), varargin{:});
lower = feval(func,out(2), varargin{:});
while((out(3)-out(2)) > delta)
    if (upper > lower)
        out(4) = out(3);
        out(3) = out(2);
        width = out(4)-out(1);
        out(2) = out(4)-gr*width;
        upper = lower;
        lower = feval(func,out(2), varargin{:});
    else
        out(1) = out(2);
        out(2) = out(3);
        width = out(4)-out(1);
        out(3) = out(1)+width*gr;
        lower = upper;
        upper = feval(func,out(3), varargin{:});
    end
end
min = out(2);
minval = feval(func,out(2), varargin{:});

