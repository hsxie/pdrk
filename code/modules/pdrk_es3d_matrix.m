% 18-10-06 08:23 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO),
% etc ...
% This file set the pdrk kernel matrix elements
% 18-10-13 10:37 update to with loss-cone distribution
% This file is for electrostatic 3D case. 2018-10-20 18:14

% % -- main program begin to set the matrix elements --
k=sqrt(kz^2+kx^2);
bsab=kx*rhocsab;
bsab(abs(bsab)<1e-50)=1e-50;  % to avoid singular when k_perp=0

bsab2=bsab.^2;
M=sparse(NN,NN);
snj=0;

% initialize
csnj=zeros(1,3*SNJ);
bsnj=csnj.*0;

for s=1:S % species
  for n=-N:N % Bessel function
    %             Gamn=exp(-bs2(s))*besseli(n,bs2(s)); % large k_perp will NaN
    %             Gamnp=exp(-bs2(s))*(besseli(n+1,bs2(s))+...
    %                 besseli(n-1,bs2(s))-2*besseli(n,bs2(s)))/2;
    
    for j=1:J % poles of Z(zeta)
      snj=snj+1;
      
      for iab=1:2 % 2018-10-13 10:58
      
      Gamn=besseli(n,bsab2(iab,s),1); % 2014-10-13 12:51
      Gamnp=(besseli(n+1,bsab2(iab,s),1)+...
        besseli(n-1,bsab2(iab,s),1)-2*besseli(n,bsab2(iab,s),1))/2;
  
      csnj(snj)=czj(j)*kz*vtzs(s)+kz*vds(s)+n*wcs(s); %

      cnj=csnj(snj);
      
%       bsnj(snj)=bsnj(snj)+rsab(iab)*Gamn*kDs(s)^2/k^2*bzj(j)*(...
%           lmdTab(iab,s)*n*wcs(s)+czj(j)*kz*vtzs(s)); % wrong rsab(iab)
      bsnj(snj)=bsnj(snj)+rsab(iab,s)*Gamn*kDs(s)^2/k^2*bzj(j)*(...
          lmdTab(iab,s)*n*wcs(s)+czj(j)*kz*vtzs(s)); % should be rsab(iab,s), 18-10-21 17:06
 
      end
    end
  end
end

for snj=1:SNJ % set the eigen matrix
    
  % (n_snj, n_snj)
  M=M+sparse(snj,snj,csnj(snj),NN,NN);
  
  % (n_snj, E)
  M=M+sparse(snj,NN,bsnj(snj),NN,NN);

  % (E, n_snj)
  M=M+sparse(NN,snj,-csnj(snj),NN,NN);

end

% (E, E)
M=M+sparse(NN,NN,-sum(bsnj),NN,NN);

% % -- main program end of set the matrix elements --
