% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2014-08-24 00:34
% Cal_J_pole_bjcj.m, Calculate J-pole (Pade approximation) coefficient
%  bj and cj for Z_J(zeta)=sum(bj/(zeta-cj)).
% Ref: Martin1980, Ronnmark1982 p11, 
%      Xie & Xiao, PDRK Appendix A, 2014
% Solve the matrix for pj and qj, then using residue() to get bj and cj.
% This is easy, since we need just the numerical data of bj and cj. The
% analytical form of pj and qj would be more difficult to obtain, which we 
% do not care at present.
% Ref: 
%   [Xie2014] H. S. Xie & Y. Xiao, PDRK: A General Kinetic Dispersion
%             Relation Solver for Magnetized Plasma.

close all; clear; clc;

J=12; % set J-pole
I=16; % Default for Ronnmark1982 is Z_{10,6}, i.e., J=8, I=10

anm=zeros(2*J+1,1); anp=0.*anm;
for n=0:J+1                        % Z(zeta) two-side
    
    % modify anp and anm here for other distribution functions
    
     % for zeta -> 0, a_n^{+}
    anp(2*n+1)=(-1)^n*1i*sqrt(pi)/gamma(n+1); % a_{2n}
    anp(2*n+2)=-(-1)^n*gamma(1/2)/gamma(n+3/2); % a_{2n+1}
    
     % for zeta -> infty, a_n^{-}
    anm(2*n+1)=0; % a_{-2n}
    anm(2*n+2)=-gamma(n+1/2)/gamma(1/2); % a_{-(2n+1)}
end

Mpq=zeros(2*J); B=zeros(2*J,1);
% I equations for zeta -> 0; K for zeta -> infty; I+K=2J
K=2*J-I;
for j=1:I % for zeta -> 0
    B(j)=anp(j); % ap_j*q_0
    if(j<=J) % p_j
        Mpq(j,j)=1;
    end
    for k=1:j-1 % ap_k*q_{j-k}
        if(j<=J+k)
            Mpq(j,J+j-k)=-anp(k);
        end
    end
end

for j=1:K % for zeta -> infty
    B(I+j)=0; % 0
    if(j<=J) % p_j
        Mpq(I+j,J-j+1)=1; % p_{J-j}
    end
    for k=1:j % am_k*q_{J+1-j}
        if(j<2*J+k)
            Mpq(I+j,2*J-j+k)=-anm(k+1);
        end
    end
    if(j>J)
        B(I+j)=-anm(2)*1; % from q(1)=1;
    end
end
pq=Mpq\B; % solve p_j (0, 1, ..., J-1) and q_j (1, 2, ..., J)

% obtain c_j and b_j
p=pq(1:J); q=pq(J+1:2*J);
p=flipud(p); q=flipud(q); q=[q;1];
[b, c, r] = residue(p,q)

bzj=b;
czj=c;

J=length(bzj)       %
sum(bzj)            % should be '-1'
sum(bzj.*czj)       % should be '0'
sum(bzj.*czj.^2)    % should be '-1/2'

%% plot
xa=-5; xb=5; dx=(xb-xa)/100; 
ya=-3; yb=5; dy=(yb-ya)/100; 
[x,y]=ndgrid(xa:dx:xb,ya:dy:yb);
z=x+1i*y;
ZJ=0.*x;
for j=1:J
    ZJ=ZJ+bzj(j)./(z-czj(j));
end

h=figure('unit','normalized','Position',[0.01 0.37 0.5 0.5],...
    'DefaultAxesFontSize',15);

subplot(211);contour(x,y,real(ZJ)); box on;
hold on; plot(real(czj),imag(czj),'rx','Linewidth',2);
title(['Re, J=',num2str(J),', I=',num2str(I)]); 
xlabel('Re(z)'); ylabel('Im(z)');
subplot(212);contour(x,y,imag(ZJ)); box on;
hold on; plot(real(czj),imag(czj),'rx','Linewidth',2);
title(['Im, J=',num2str(J),', I=',num2str(I)]); 
xlabel('Re(z)'); ylabel('Im(z)');

print(gcf,'-dpng',['J',num2str(J),'I',num2str(I),...
    'xa',num2str(xa),'xb',num2str(xb),...
    'ya',num2str(ya),'yb',num2str(yb),'_2d.png']);

