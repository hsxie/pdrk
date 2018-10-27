% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2014-05-24 17:45
% pdrk_es3d.m, Plasma Dispersion Relastion solver, Kinetic -> ES magnetized 
% -> Harris dispersion relation (Ref: Gurnett, 2005, p387).
% bi-Maxwellian equilibrium distribution.
% transform to matrix eigenvalue problem lambda*X=M*X
% J-pole approximation for Z(zeta)=sum(b_j/(zeta-c_j))
% Ref: 
%   [Xie2014] H. S. Xie & Y. Xiao, PDRK: A General Kinetic Dispersion
%             Relation Solver for Magnetized Plasma.

close all; clear; clc;

N=10; B0=1.0;

%% read input parameters
par = importdata('pdrk_es3d.in', ' ', 1); 
[S, col]=size(par.data);
if(col~=6)
    disp('Wrong input data !!!');
end

for s=1:S
    qs(s)=par.data(s,1); % charge
    ms(s)=par.data(s,2); % mass
    ns0(s)=par.data(s,3); % desity
    Tzs(s)=par.data(s,4); % parallel temperature
    Tps(s)=par.data(s,5); % perp temperature
    vs0(s)=par.data(s,6); % para drift velocity
end
vtzs=sqrt(2*Tzs./ms); % para thermal velocity
lambdaDs=sqrt(Tzs./(ns0.*qs.^2)); % Debye length, Tzs
kDs=1./lambdaDs;
wps=sqrt(ns0.*qs.^2./ms); % plasma frequency
wcs=B0*qs./ms; % cyclotron frequency
rhocs=sqrt(Tps./ms)./abs(wcs); % cyclotron radius

J=8; % J-pole
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
elseif(J==12) % from Cal_J_pole_bjcj.m
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
J=length(bzj)
sum(bzj)
sum(bzj.*czj)
sum(bzj.*czj.^2)


SNJ=S*(2*N+1)*J;
kk=[]; ww=[]; kxx=[]; kzz=[];
tic;
for kj=0.05:0.1:5
    kz=0*kj;
    kx=kj;
    
    k=sqrt(kz^2+kx^2);
    
    bs=kx*rhocs;
    bs2=bs.^2;
    M=zeros(SNJ,SNJ);
    snj=0;
    aSNJ=1;
    for s=1:S % species 
        aSNJ=aSNJ+(kDs(s)^2/k^2);
        for n=-N:N % Bessel function
%             Gamn=exp(-bs2(s))*besseli(n,bs2(s)); % large k_perp will NaN
            Gamn=besseli(n,bs2(s),1); % 2014-10-13 12:51

            for j=1:J  % poles of Z(zeta)
                snj=snj+1;
                tmp=n*wcs(s)+czj(j)*kz*vtzs(s);
                csnj(snj)=kz*vs0(s)+tmp;
                bsnj(snj)=(kDs(s)^2/k^2)*Gamn*bzj(j)*(tmp-...
                    (1-Tzs(s)/Tps(s))*n*wcs(s));
                aSNJ=aSNJ+(kDs(s)^2/k^2)*Gamn*(bzj(j));
            end
        end
    end

    for snj=1:SNJ % set the eigen matrix
        M(snj,:)=-bsnj(snj);
        M(snj,snj)=-bsnj(snj)+aSNJ*csnj(snj);
    end

    d=eig(M); % solve the roots

    omega=d;
    [wi,ind]=sort(imag(omega),'descend');
%     [wi,ind]=sort(abs(real(omega)),'ascend');
    w=omega(ind);
    
    kk=[kk,k]; kxx=[kxx,kx]; kzz=[kzz,kz];
    ww=[ww,w];
end
runtime=toc;

%%
h=figure('unit','normalized','Position',[0.01 0.4 0.65 0.4],...
    'DefaultAxesFontSize',15);
kkp=kxx;
subplot(121);
for n=0:6
    plot([0,max(kkp)],[n*abs(wcs(1)),n*abs(wcs(1))],'g--','LineWidth',2);
end
for j=1:SNJ
%     ind=find(imag(ww)<-1e-3); ww(ind)=NaN;   
      % change here to remove artificial solutions
%     if(min((imag(ww(j,:))))>-4e-0) 
        subplot(121);
        plot(kkp,real(ww(j,:)),'r*','LineWidth',2);
        hold on;
        subplot(122);
        plot(kkp,imag(ww(j,:)),'r*','LineWidth',2);
        hold on;
%     end
end

subplot(121);
title(['\omega_p/\omega_c=',num2str(abs(wps(1)/wcs(1))),...
    ', N=',num2str(N)]);
xlim([0,max(kkp)]); 
ylim([0,5.5]);
xlabel('k_{\perp}\rho_c'); ylabel('\omega_r/\omega_c');
subplot(122);
xlabel('k_{\perp}\rho_c'); ylabel('\omega_i/\omega_c');
xlim([0,max(kkp)]);

print('-dpng',['pdrk_es3d_wp',num2str(abs(wps(1)/wcs(1))),...
    '_N',num2str(N),'_J',num2str(J),'.png']);


