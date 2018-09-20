%Erin Carson
%lejapointscomplex.m
%Edited 5/13/2015

%Input:
%a: vector containing real parts of eigenvalues
%b: vector containing imaginary parts of eigenvalues. So a(j) + i*b(j) are
%eigenvalues of A
%s: number of leja points returned

%Output:
%ljapnts: vector of s Leja points to be used in basis computation (complex)

function ljapnts = lejapointscomplex(a,b,s)

%Set number of points per "spoke" (see Philippe and Reichel, 2012)
n_perspoke = 10;

i = sqrt(-1);

%Find center of eigenvalues
a_centr = mean(a);
b_centr = mean(b);

%For each eigenvalue, select points along the "spoke" line, running from
%the eigenvalue to the center of all eigenvalues
for j = 1:numel(a)
     
    dist_y = b(j) - b_centr;
    y_step = dist_y/n_perspoke;
    y_pnts = zeros(1,n_perspoke);
    y_pnts(1) = b_centr + y_step;
    for i = 2:numel(y_pnts)
        y_pnts(i) = y_pnts(i-1)+y_step;
    end
        
    dist_x = a(j) - a_centr;
    x_step = dist_x/n_perspoke;
    x_pnts = zeros(1,n_perspoke);
    x_pnts(1) = a_centr + x_step;
    for i = 2:numel(x_pnts)
        x_pnts(i) = x_pnts(i-1)+x_step;
    end
    
    if(numel(y_pnts) == 0)
        y_pnts = zeros(1,n_perspoke);
    end
    if(numel(x_pnts) == 0)
        x_pnts = zeros(1,n_perspoke);
    end
    y_pnts = y_pnts(1:n_perspoke);
    x_pnts = x_pnts(1:n_perspoke);
    pnts((n_perspoke*(j-1)+1):(n_perspoke*(j-1)+n_perspoke)) = x_pnts + i*y_pnts;
       
end

pnts = unique(pnts);
yy = lejaorder(pnts,s);

ljapnts = yy;

end

function z = lejaorder(x,n)

[v, ind ] = max(abs((x)));
f_ind = ind;%round(numel(x)*rand(1));
z(1) = x(f_ind);

k = 1;
if( real(z(1)) ~= z(1) )
    if(imag(z(1)) > 0)
        z(2) = conj(z(1));
    else
        
        z(1) = conj(z(1));
        z(2) = conj(z(1));
    end
    inds = find(x==z(2));
    x(inds) = [];
    inds = find(x==z(1));
    x(inds) = [];
    k = 2;
end

tmp = abs(x-z(1));




while k <= n
   for ii = 1:numel(x)
        tmp(ii) = tmp(ii)*abs(x(ii)-z(k));
   end
   
    [v, index] = max(abs(tmp));
  
    %If we haven't used this point yet
    if(numel(find(z==x(index))) == 0)
        %Find if guy is anywhere else in the list, remove those
        %entries
        z(k+1) = x(index);
        x(index) = [];
        tmp(index) = [];
        
   
        %If there is an imaginary component
        if( real(z(k+1)) ~= z(k+1) )
            
            %Find if conjugate is anywhere else in the list, remove those
            %entries
            inds = find(x==conj(z(k+1)));
            x(inds) = [];
            tmp(inds) = [];
           
            %If the imaginary component is positive
            %put its conjugate before it
            
            if(imag(z(k+1))>0)               
                z(k+2) = conj(z(k+1));
            else
                z(k+1) = conj(z(k+1));
                z(k+2) = conj(z(k+1));
            end

            %Update list for next search
            tmp = tmp.*abs(x-(ones(1,numel(x)).*z(k+1)));
            
            %
            k  = k + 1;
            
        end
        k = k+1;
    end
    
end

z = z(1:n);
if(abs(imag(z(end))) > 1e-10 && z(end) ~= conj(z(end-1)))
    z(end) = real(z(end));
end

end

