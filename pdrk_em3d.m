% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2014-06-01 17:00
% pdrk_em3d.m, Plasma Dispersion Relastion solver, kinetic, EM3D,
% bi-Maxwellian equilibrium distribution with parallel drift.
% Transform to matrix eigenvalue problem lambda*X=M*X.
% J-pole approximation for Z(zeta)=sum(b_j/(zeta-c_j))
%
% sparse, eigs
%
% 2014-08-22 01:09, support T_perp/T_para \neq 1,  vs0 \neq 0
%
% Ref:
%  [Xie2016] H. S. Xie & Y. Xiao, PDRK: A General Kinetic Dispersion
%    Relation Solver for Magnetized Plasma, Plasma Science and Technology,
%    Vol.18, No.2, p97 (2016). (Also 2014arXiv, note: have several typos)
%
% Documents/codes/Erratum: http://hsxie.me/codes/pdrk/
%
% 18-09-28 12:42 We find that omega_cs=-qB/m in Miyamoto2004, not the
% standard omega_cs=qB/m. This sign difference affect several terms. This
% version we have fixed the bugs of the sign of omega_cs.
% 18-10-01 09:43 based on new derivations, fixed a bug in b33 term, where a
% wcs^2 is missed.
% 16:01 Benchmark with Gary1993 Fig.7.4 and several other cases ok.
%
% Ref (typos/bugs fixed):
%  H. S. Xie, Detailed Derivations of PDRK-EM3D Equations, 2018-10-03 (
%  10 pages).
% 
% 18-10-03 13:14 This version should be bugs free now.

close all;
clear; clc;

% enlarge N to make sure the results are convergent
N=2; % Number harmonics
J=8; % J-pole, usually J=8 is sufficient

B0=100.0E-9; c2=(2.9979E8)^2; epsilon0=8.8542E-12;
mu0=1/(c2*epsilon0);
kB=1.3807e-23;

% % read input parameters
par = importdata('pdrk_em3d.in', ' ', 1);
[S, col]=size(par.data);
if(col~=6)
  disp('Wrong input data !!!');
end

for s=1:S
  qs0(s)=par.data(s,1); % charge
  ms0(s)=par.data(s,2); % mass
  ns0(s)=par.data(s,3); % desity, m^-3
  Tzs0(s)=par.data(s,4); % parallel temperature, eV
  Tps0(s)=par.data(s,5); % perp temperature, eV
  vs0(s)=par.data(s,6); % para drift velocity
end

qs=qs0*1.6022e-19; % * electron charge, e -> C (coulomb)
ms=ms0*9.1094E-31; % * electron mass, m_e -> kg
Tzs=Tzs0*1.6022e-19/kB; % T//, eV -> K (eV -> J * J -> K)
Tps=Tps0*1.6022e-19/kB; % Tpr, eV -> K
vs0=vs0*sqrt(c2); % vs0, speed of light c -> m/s

vtzs=sqrt(2*kB*Tzs./ms); % para thermal velocity, note the sqrt(2)
lambdaDs=sqrt(epsilon0*kB*Tzs./(ns0.*qs.^2)); % Debye length, Tzs
kDs=1./lambdaDs;
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wcs=B0*qs./ms; % cyclotron frequency
% rhocs=sqrt(kB*Tps./ms)./abs(wcs); % cyclotron radius
rhocs=sqrt(kB*Tps./ms)./wcs; % cyclotron radius, 2018-06-13 21:47

wps2=wps.^2;
lmdT=Tzs./Tps;

betasz=2*mu0*kB.*ns0.*Tzs./B0^2; % beta_para
betasp=2*mu0*kB.*ns0.*Tps./B0^2; % beta_perp
vA=B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed
cS=sqrt(2*min(kB*Tzs)/max(ms));

% change here if the second species is not ion
wci=abs(wcs(2)); % omega_{ci}
wpi=sqrt(ns0(1)*qs(2)^2/ms(2)/epsilon0); % ion plasma frequency using electron density
cwp=sqrt(c2)/wpi; % c/omega_{pi}
vAwp=vA/wci; % v_A/omega_{ci}

if(J==8)
  % Ronnmark1982, 8-pole for Z function, and Z'
  bzj(1)=-1.734012457471826E-2-4.630639291680322E-2i;
  bzj(2)=-7.399169923225014E-1+8.395179978099844E-1i;
  bzj(3)=5.840628642184073+9.536009057643667E-1i;
  bzj(4)=-5.583371525286853-1.120854319126599E1i;
  czj(1)=2.237687789201900-1.625940856173727i;
  czj(2)=1.465234126106004-1.789620129162444i;
  czj(3)=.8392539817232638-1.891995045765206i;
  czj(4)=.2739362226285564-1.941786875844713i;
  
  bzj(5:8)=conj(bzj(1:4));
  czj(5:8)=-conj(czj(1:4));
elseif(J==12) % from Cal_J_pole_bjcj.m, Xie2016
  bzj(1)=    -0.00454786121654587 - 0.000621096230229454i;
  bzj(2)=    0.215155729087593 + 0.201505401672306i;
  bzj(3)=    0.439545042119629 + 4.16108468348292i;
  bzj(4)=    -20.2169673323552 - 12.8855035482440i;
  bzj(5)=    67.0814882450356 + 20.8463458499504i;
  bzj(6)=    -48.0146738250076 + 107.275614092570i;
  
  czj(1)=    -2.97842916245164 - 2.04969666644050i;
  czj(2)=    2.25678378396682 - 2.20861841189542i;
  czj(3)=    -1.67379985617161 - 2.32408519416336i;
  czj(4)=    -1.15903203380422 - 2.40673940954718i;
  czj(5)=    0.682287636603418 - 2.46036501461004i;
  czj(6)=    -0.225365375071350 - 2.48677941704753i;
  
  bzj(7:12)=conj(bzj(1:6));
  czj(7:12)=-conj(czj(1:6));
elseif(J==4) % Martin1980
  bzj(1)=0.546796859834032 + 0.037196505239277i;
  bzj(2)=-1.046796859834027 + 2.101852568038518i;
  czj(1)=1.23588765343592 - 1.21498213255731i;
  czj(2)=-0.378611612386277 - 1.350943585432730i;
  bzj(3:4)=conj(bzj(1:2));
  czj(3:4)=-conj(czj(1:2));
end
% J=length(bzj);
% sum(bzj) % should be -1
% sum(bzj.*czj) % should be 0
% sum(bzj.*czj.^2) % should be -1/2

SNJ=S*(2*N+1)*J;
SNJ1=SNJ+1;
SNJ3=3*SNJ1;
NN=SNJ3+6;

% modified here if you only need one or several solutions around initial
% guess value wg
sp=0; % sparse eigs() or eig()
if(sp==0)
    nw=NN; % number of roots to obtain
else % using sparse matrix
    nw=1; % only nw solutions around the initial guesses are given
    wg=1.5; % inital guess
end

if(sp==1)
  wg=wg0*wci;
end

kk=[]; kxx=[]; kzz=[]; ww=[];
tic;
% modified here if you want to scan other parameters instead of k
for kj=0.01:0.02:1.0 % kj=0.0:2e-6:2.5e-4 %
    
    theta=71*pi/180; % the angle between k and B0
    kz=cos(theta)*kj/cwp; kx=sin(theta)*kj/cwp;
    
    % % -- main program begin to set the matrix elements --
    k=sqrt(kz^2+kx^2);
    bs=kx*rhocs;
    bs(abs(bs)<1e-50)=1e-50;  % to avoid singular when k_perp=0

    bs2=bs.^2;
    M=sparse(NN,NN);
    snj=0;

    % initialize
    b11=0; b12=0; b13=0; b21=0; b22=0; b23=0; b31=0; b32=0; b33=0;
    csnj=zeros(1,3*SNJ);
    b11snj=csnj.*0; b12snj=csnj.*0; b13snj=csnj.*0;
    b21snj=csnj.*0; b22snj=csnj.*0; b23snj=csnj.*0;
    b31snj=csnj.*0; b32snj=csnj.*0; b33snj=csnj.*0;

    for s=1:S % species
      for n=-N:N % Bessel function
        %             Gamn=exp(-bs2(s))*besseli(n,bs2(s)); % large k_perp will NaN
        %             Gamnp=exp(-bs2(s))*(besseli(n+1,bs2(s))+...
        %                 besseli(n-1,bs2(s))-2*besseli(n,bs2(s)))/2;
        Gamn=besseli(n,bs2(s),1); % 2014-10-13 12:51
        Gamnp=(besseli(n+1,bs2(s),1)+...
          besseli(n-1,bs2(s),1)-2*besseli(n,bs2(s),1))/2;
        for j=1:J % poles of Z(zeta)
          snj=snj+1;

          csnj(snj)=czj(j)*kz*vtzs(s)+kz*vs0(s)+n*wcs(s);

          cnj=csnj(snj);
          bj0=vs0(s)+(1-1/lmdT(s))*czj(j)*vtzs(s);

          %
          if(n==0)  % for A_nj
            bnj1=bj0/(czj(j)*vtzs(s)+vs0(s)); % avoid cnj=0
          else
            bnj1=kz*bj0/cnj;
          end
          bnj2=1-bnj1;

          tmp=wps2(s)*bzj(j);

          b11snj(snj)=tmp*bnj2*n^2*Gamn/bs2(s);
          b11=b11+tmp*bnj1*n^2*Gamn/bs2(s);

          b12snj(snj)=tmp*bnj2*1i*n*Gamnp;
          b12=b12+tmp*bnj1*1i*n*Gamnp;
          b21snj(snj)=-b12snj(snj);
          b21=-b12;

          b22snj(snj)=tmp*bnj2*(n^2*Gamn/bs2(s)-2*bs2(s)*Gamnp);
          b22=b22+tmp*bnj1*(n^2*Gamn/bs2(s)-2*bs2(s)*Gamnp);

          %
          if(n==0)  % for eta_n*A_nj
            bnj1=0; % avoid cnj=0 when kz=0
          else
            bnj1=n*wcs(s)*bj0/cnj/vtzs(s);
          end
          bnj2=czj(j)/lmdT(s)+bnj1;

          b13snj(snj)=tmp*bnj2*n*sqrt(2*lmdT(s))*Gamn/bs(s);
          b13=b13-tmp*bnj1*n*sqrt(2*lmdT(s))*Gamn/bs(s);
          b31snj(snj)=b13snj(snj);
          b31=b13;

          b23snj(snj)=-1i*tmp*bnj2*sqrt(2*lmdT(s))*Gamnp*bs(s);
          b23=b23+1i*tmp*bnj1*sqrt(2*lmdT(s))*Gamnp*bs(s);
          b32snj(snj)=-b23snj(snj);
          b32=-b23;

          %
          if(bj0==0 || kz==0)  % for eta_n^2*A_nj
            bnj1=0;
            bnj2=czj(j)*czj(j);
          else
            % bnj1=n^2*bj0/cnj/vtzs(s)^2/kz;
            bnj1=n^2*wcs(s)^2*bj0/cnj/vtzs(s)^2/kz; % !!fixed bug of missed wcs^2, 18-10-01 10:19
            bnj2=(vs0(s)/vtzs(s)+czj(j))*czj(j)/lmdT(s)+...
              n*wcs(s)*bj0*(1-n*wcs(s)/cnj)/vtzs(s)^2/kz;
          end

          b33snj(snj)=tmp*bnj2*2*lmdT(s)*Gamn;
          b33=b33+tmp*bnj1*2*lmdT(s)*Gamn;
        end
      end
    end

    for snj=1:SNJ % set the eigen matrix
      jjx=snj+0*SNJ1;
      jjy=snj+1*SNJ1;
      jjz=snj+2*SNJ1;
      % v_snjx
      M=M+sparse(jjx,jjx,csnj(snj),NN,NN)+...
        sparse(jjx,SNJ3+1,b11snj(snj),NN,NN)+...
        sparse(jjx,SNJ3+2,b12snj(snj),NN,NN)+...
        sparse(jjx,SNJ3+3,b13snj(snj),NN,NN);

      % v_snjy
      M=M+sparse(jjy,jjy,csnj(snj),NN,NN)+...
        sparse(jjy,SNJ3+1,b21snj(snj),NN,NN)+...
        sparse(jjy,SNJ3+2,b22snj(snj),NN,NN)+...
        sparse(jjy,SNJ3+3,b23snj(snj),NN,NN);

      % v_snjz
      M=M+sparse(jjz,jjz,csnj(snj),NN,NN)+...
        sparse(jjz,SNJ3+1,b31snj(snj),NN,NN)+...
        sparse(jjz,SNJ3+2,b32snj(snj),NN,NN)+...
        sparse(jjz,SNJ3+3,b33snj(snj),NN,NN);

    end

    % E(J), J_{x,y,z}=j_{x,y,z}+sum(v_snj{x,y,z})
    tp=-1;
    jj=(0*SNJ1+1):(1*SNJ1); ii=0.*jj+SNJ3+1; M=M+sparse(ii,jj,tp,NN,NN);
    jj=(1*SNJ1+1):(2*SNJ1); ii=0.*jj+SNJ3+2; M=M+sparse(ii,jj,tp,NN,NN);
    jj=(2*SNJ1+1):(3*SNJ1); ii=0.*jj+SNJ3+3; M=M+sparse(ii,jj,tp,NN,NN);

    % jx(E), jy(E), jz(E)
    M=M+sparse(1*SNJ1,SNJ3+1,b11,NN,NN)+...
      sparse(1*SNJ1,SNJ3+2,b12,NN,NN)+...
      sparse(1*SNJ1,SNJ3+3,b13,NN,NN)+...
      sparse(2*SNJ1,SNJ3+1,b21,NN,NN)+...
      sparse(2*SNJ1,SNJ3+2,b22,NN,NN)+...
      sparse(2*SNJ1,SNJ3+3,b23,NN,NN)+...
      sparse(3*SNJ1,SNJ3+1,b31,NN,NN)+...
      sparse(3*SNJ1,SNJ3+2,b32,NN,NN)+...
      sparse(3*SNJ1,SNJ3+3,b33,NN,NN);

    % E(B)
    M=M+sparse(SNJ3+1,SNJ3+5,c2*kz,NN,NN)+...
      sparse(SNJ3+2,SNJ3+4,-c2*kz,NN,NN)+...
      sparse(SNJ3+2,SNJ3+6,c2*kx,NN,NN)+...
      sparse(SNJ3+3,SNJ3+5,-c2*kx,NN,NN);

    % B(E)
    M=M+sparse(SNJ3+4,SNJ3+2,-kz,NN,NN)+...
      sparse(SNJ3+5,SNJ3+1,kz,NN,NN)+...
      sparse(SNJ3+5,SNJ3+3,-kx,NN,NN)+...
      sparse(SNJ3+6,SNJ3+2,kx,NN,NN);
  
    % % -- main program end of set the matrix elements --
    
    % calculate the solution use either eig() or eigs()
    if(sp==0) % solve all the roots, slow for large N
      % d=eig(full(M));
      % it's found in pdrf (fluid version) that vpa(eig) may more accurate
      d0=vpa(eig(full(M)),16);d=double(d0);
    else % solve only nw solutions, fast for large N, similar to WHAMP
      % d=eigs(M,nw,wg);
      d0=vpa(eigs(M,nw,wg),16);d=double(d0);
    end
    omega=d;
    [wi,ind]=sort(imag(omega),'descend'); % sort the solution by growth rate
    w=omega(ind);
    wg=w(1); % update the new initial guess use the previous solution

    kk=[kk,k];
    kxx=[kxx,kx]; kzz=[kzz,kz];
    ww=[ww,w];

end %jkn
runtime=toc;

%% modify here to plot the solutions
close all;
kkp=kk;

h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.3],...
    'DefaultAxesFontSize',15);
subplot(121);

for n=1:1
%     plot([0,max(kk)],[n*abs(wci),n*abs(wci)],'c--','LineWidth',2);
    hold on;
end

kkkp=reshape(repmat(kkp,nw,1),1,[]);
www=reshape(ww,1,[]);
nkk=length(kkkp);
jj=[];
for j=1:nkk % to filter some extraneous solutions
%     if(imag(www(j))<-2.0e-0 || imag(www(j))>1.0e-0)
%     if(min(imag(www(j,:)))/max(abs(real(www(j,:))))>-0.96e-0 ...
%             && max(abs(real(www(j,:))))/max(kk)<2.902/10  || min(imag(www(j,:)))>-7.5e-1)
%     if(abs(imag(www(j)))/kkkp(j)>2.0e-1 || abs(imag(www(j))/kkkp(j)+0.17)<1.0e-3) 
%         jj=[jj,j];
%     end
end
www(jj)=[];
kkkp(jj)=[];

plt=1; %
if(plt==1)
    for j=1:1 % plot the first several largest growth rate modes
    %for j=1:5
    subplot(121);
%     plot(kkkp,real(www/wci),'bx','LineWidth',2); hold on;
    plot(kk(1,:)*cwp,real(ww(j,:)/wci),'bx','LineWidth',2); hold on;
    xlabel('kc/\omega_{p}'); ylabel('\omega_r/\omega_{cp}'); box on;
    title(['(a) \beta_{||}=',num2str(betasz,3),', \beta_\perp=',num2str(betasp,3)]);
    ylim([-0.001,0.001]);
    %ylim([-5.0,25.0]);
    subplot(122);
%     plot(kkkp,imag(www/wci),'bx','LineWidth',2); hold on;
    plot(kk(1,:)*cwp,imag(ww(j,:)/wci),'bx-','LineWidth',2); hold on;
    xlabel('kc/\omega_{p}'); ylabel('\omega_i/\omega_{cp}'); box on;
    title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', \theta=',...
        num2str(theta*180/pi),'^{\circ}, (N=',num2str(N),', J=',num2str(J),')']);
    % v_{s0}/v_A=',num2str(vs0/vA,3)
%     ylim([-0.10,0.1]);
    end
else
    ind=find(abs(imag(www)/wci)<1e-5);
    ind1=find(imag(www)/wci<-0.01); % to filter some extraneous solutions
    www1=www; www1(ind)=NaN+1i*NaN;www1(ind1)=NaN+1i*NaN;
    subplot(121);
    plot(kkkp*cwp,real(www1)/wci,'g+','LineWidth',2); hold on;
    xlabel('kc/\omega_{cp}'); ylabel('\omega_r/\omega_{cp}'); box on;
    title(['(a) \beta_{||}=',num2str(betasz,3),', \beta_\perp=',num2str(betasp,3)]);
    %ylim([-0.7,0.7]);
    ylim([-0.25,0.25]);
%     ylim([0,100]);
    subplot(122);
    plot(kkkp*cwp,imag(www1)/wci,'g+','LineWidth',2); hold on;
    xlabel('kc/\omega_{cp}'); ylabel('\omega_i/\omega_{cp}'); box on;
    title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', \theta=',...
        num2str(theta*180/pi),'^{\circ}, (N=',num2str(N),', J=',num2str(J),')']);
    % v_{s0}/v_A=',num2str(vs0/vA,3)
    ylim([-0.01,0.01]);
end

% 2018-04-14, R. Denton's WHAMP result of Gary1993 book Fig.7.4, benchmark
k_whamp=0:.05:.60;
gamma_whamp=[0 1.797E-03 3.525E-03 5.107E-03 6.481E-03  7.55E-03 ...
     8.29E-03  8.55E-03 8.277E-03 7.368E-03    5.720E-03    3.230E-03 0];
plot(k_whamp,gamma_whamp,'r+','linewidth',2);

legend('pdrk-xie','whamp-denton','Location','best'); legend('boxoff');

%%
% figure;
% 
% plot3(real(www),imag(www),kkkp,'r+','LineWidth',2); hold on;
% xlabel('Re'); ylabel('Im'); zlabel('k');

%%

% %%
% fid = fopen('w2.dat','wt');
% for i=1:length(kplot)
%   fprintf(fid,'%10g\t',kplot(i));
%   fprintf(fid,'%15g\t',real(wplot(i)));
%   fprintf(fid,'%10g\n',imag(wplot(i)));
% end
% fclose(fid);

%%

print('-dpng',['pdrk_em3d_wp',num2str(abs(wps(1)/wcs(1))),...
  '_N',num2str(N),'_J',num2str(J),'_sp',num2str(sp),...
  '_nw',num2str(nw),'_vd',num2str(max(abs(vs0))),'_lmdT',...
  num2str(max(abs(lmdT))),'_theta=',num2str(theta*180/pi),'_new.png']);
