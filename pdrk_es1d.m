% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-04-18 19:09
% 2014-05-24 22:40, rewrite
% pdrk_es1d.m, Plasma Dispersion Relastion solver, Kinetic, ES1D
% transform to matrix eigenvalue problem lambda*X=MX
% Landau damping: 1+sum_s{[1+zeta_s*Z(zeta_s)]/(k*lamda_Ds)^2}=0,
% zeta_s=(w-k*vs0)/(k*vts)
% J-pole approximation for Z_J(zeta)=sum(bj/(zeta-cj))
% Ref: 
%   [Xie2014] H. S. Xie & Y. Xiao, PDRK: A General Kinetic Dispersion
%             Relation Solver for Magnetized Plasma.

close all; clear; clc;

%% read input parameters
par = importdata('pdrk_es1d.in', ' ', 1); 
[S, col]=size(par.data);
if(col~=5)
    disp('Wrong input data !!!');
end

for j=1:S
    qs(j)=par.data(j,1);
    ms(j)=par.data(j,2);
    ns0(j)=par.data(j,3);
    Ts(j)=par.data(j,4);
    vs0(j)=par.data(j,5);
end
vts=sqrt(2*Ts./ms); % thermal velocity
lambdaDs=sqrt(Ts./(ns0.*qs.^2)); % Deybe length
kDs=1./lambdaDs;
wps=sqrt(ns0.*qs.^2./ms); % plasma frequency

J=8;
if(J==8)
    opt=1;
    if(opt==1)
        % Ronnmark1982, 8-pole for Z function, and Z'
        bzj(1)=-1.734012457471826E-2-4.630639291680322E-2i;
        bzj(2)=-7.399169923225014E-1+8.395179978099844E-1i;
        bzj(3)=5.840628642184073+9.536009057643667E-1i;
        bzj(4)=-5.583371525286853-1.120854319126599E1i;
        czj(1)=2.237687789201900-1.625940856173727i;
        czj(2)=1.465234126106004-1.789620129162444i;
        czj(3)=.8392539817232638-1.891995045765206i;
        czj(4)=.2739362226285564-1.941786875844713i;
    else % new calculation
        bzj(1)=  -0.0173401116032742 - 0.0463064419344598i;
        bzj(2)=  -0.739917851897683 + 0.839518298070637i;
        bzj(3)=  5.84063227513760 + 0.953602843950785i;
        bzj(4)=  -5.58337431170864 - 11.2085508179677i;
        czj(1)=   2.23768772215616 - 1.62594103256666i;
        czj(2)=   1.46523409042510 - 1.78962030806222i;
        czj(3)=   0.839253965702731 - 1.89199521968963i;
        czj(4)=   0.273936217871668 - 1.94178704551807i;
    end

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
elseif(J==4)
    opt=2;
    if(opt==1)
        % Martin1980, 4-pole
        bzj(1)=0.5468-0.0372i;
        bzj(2)=-1.0468+2.1018i;
        czj(1)=-1.2359-1.2150i;
        czj(2)=-0.3786-1.3509i;
    else % new calculation
        bzj(1)=0.546796859834032 + 0.037196505239277i;
        bzj(2)=-1.046796859834027 + 2.101852568038518i;
        czj(1)=1.23588765343592 - 1.21498213255731i;
        czj(2)=-0.378611612386277 - 1.350943585432730i;
    end
    
    bzj(3:4)=conj(bzj(1:2));
    czj(3:4)=-conj(czj(1:2));
elseif(J==3)
    % Martin1980, 3-pole
    bzj(1)=0.1822+0.5756i;
    bzj(2)=-1.3643;
    czj(1)=-0.9217-0.9091i;
    czj(2)=-1.0204i;
    
    bzj(3)=conj(bzj(1));
    czj(3)=-conj(czj(1));
    
elseif(J==2)
    % Huba2009, 2-pole
    bzj(1)=-(0.5+0.81i);
    czj(1)=0.51-0.81i;
    bzj(2)=conj(bzj(1));
    czj(2)=-conj(czj(1));   
    
else
    % Martin1979, 2-pole
    bzj(1)=-(0.5+1.2891i);
    czj(1)=0.5138-1.0324i;
    bzj(2)=conj(bzj(1));
    czj(2)=-conj(czj(1));
end
J=length(bzj)
sum(bzj)
sum(bzj.*czj)
sum(bzj.*czj.^2)


%%
kk=[]; ww=[];
for k=0.4:-0.02:0.01

    SJ=S*J;

    sj=0;
    M=zeros(SJ,SJ);
    for s=1:S
        for j=1:J;
            sj=sj+1;
            csj(sj)=k*(czj(j)*vts(s)+vs0(s));
            bsj(sj)=vts(s)*bzj(j)*czj(j)*(kDs(s)^2/k);
        end
    end

    for sj=1:SJ
        M(sj,:)=-bsj(sj);
        M(sj,sj)=-bsj(sj)+csj(sj);
    end

    d=eig(M);
    omega=d;
    [wi,ind]=sort(imag(omega),'descend');
    w=omega(ind);
    
    kk=[kk,k];
    ww=[ww,w];
end

%%
h=figure('unit','normalized','Position',[0.01 0.57 0.6 0.35],...
    'DefaultAxesFontSize',15);

subplot(121);
plot(kk,imag(ww(1,:)),'b-',kk,imag(ww(2,:)),'g--',kk,imag(ww(3,:)),...
    'c--','LineWidth',2);
xlabel('k'); ylabel('\omega_i'); 
% title(['(a) growth rate, n_b=',num2str(ns0(2))]); %
ylim([-0.5,0.3]);

subplot(122);
plot(kk,real(ww(1,:)),'b--',kk,real(ww(2,:)),'g.',kk,real(ww(3,:)),...
    'c.','LineWidth',2);
% title(['(b) frequency, T_b=T_e, v_{b}=',num2str(vs0(2)/vts(1)),'v_{te}']);
xlabel('k'); ylabel('\omega_r'); 
ylim([-1.5,3]);

% for jp=1:SJ % plot all solutions
%     subplot(121);
%     plot(kk,imag(ww(jp,:)),'.'); hold on;
% 
%     subplot(122);
%     plot(kk,real(ww(jp,:)),'.'); hold on;
% end

% run bumpontail;
print('-dpng',['pdrk_es1d_beam_vs_k_J',num2str(J),'.png']);
