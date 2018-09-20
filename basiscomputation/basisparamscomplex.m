% Erin Carson
% basisparamscomplex.m
% Edited 5/13/2015

%Find basis parameters to generate (s+1)-dimension Krylov basis for matrix
%A for selected basis type when A has complex eigenvalues

%Input:
%s: number of basis vectors to be generated (not counting the starting vector)
%evec: vector containing eigenvalues of A
%basis_type: string denoting which basis to use. Acceptable values are
% 'monomial', 'newton', or 'chebyshev'. If something besides these strings
% entered, will default to using monomial basis. 

%Output:
%T: an (s+1)-by-s tridiagonal matrix which stores basis parameters; represents multiplication by A in
%the chosen basis
%alp, bet, gam: vectors giving the diagonal, superdiagonal, and subdiagonals of T,
%respectively, of length s, s-1, and s, respectively. 

function [alp,bet,gam,T] = basisparamscomplex (s,evec,basis_type)

%Compute complex Leja points
ljapnts = lejapointscomplex(real(evec),imag(evec),s);

%If Newton basis requested, compute basis parameters for complex Newton
%basis (see Hoemmen, 2010)
if(strcmp(basis_type, 'newton'))
    T=zeros(s+1);
    for j = 1:s
            if imag(ljapnts(j)) > 0
                T(j,j)   = real(ljapnts(j));
                T(j,j+1) = -imag(ljapnts(j))^2;
            else
                T(j,j) = real(ljapnts(j));
            end
            T(j+1,j) = 1;
    end
    
    %Set alp, bet, gam parameters
    alp = diag(T,0);
    alp = alp(1:end-1);
    bet = diag(T,1);
    bet = bet(1:end-1);
    gam = diag(T,-1);

%If Chebyshev basis requested, compute basis parameters for complex
%Chebyshev basis
%Note: this is not the scaled and shifted version - this remains to be
%implemented
elseif(strcmp(basis_type, 'chebyshev'))
    
    if(numel(find(imag(evec)>0))==0)
        imgvals = 0;
    else
        imgvals = 1;
    end

    cc = mean(real(evec));
    aa=abs(max(abs(evec))-cc);
    bb=max(abs(imag(evec)));
    dd=sqrt(aa^2-bb^2);
    gamm = max(aa,bb);
          
    %Set alp, bet, gam parameters
    alp = cc.*ones(s,1);
    bet = (dd)^2/(4*gamm).*ones(s-1,1);
    gam = [2*gamm; gamm.*ones(s-1,1)];
    
else %monomial basis; note - should not reach this point if basis_type is set as a valid option; this is only a fallback
    alp = zeros(s,1);
    bet = zeros(s-1,1);
    gam = ones(s,1);
end

%Set diagonals of T matrix
T = zeros(s+1,s);
T = diag(alp,0) + diag(bet,1) + diag(gam(1:end-1),-1);
T=[T;zeros(1,s-1),gam(end)];

