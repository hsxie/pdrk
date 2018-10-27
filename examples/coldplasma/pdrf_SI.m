% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-05-03 22:02
% pdrf.m (Plasma Dispersion Relation - Fluid version), general dispersion
% relation solver for multi-fluid plasma, ver 1.0. Full-matrix method.
%
% B0=(0,0,B0), k=(kx,0,kz), vj0=(vj0x,vj0y,vj0z),
% dln(nj0)=(epsnjx,epsnjy,0)
%
% Ref: 
%  [Xie2014] H. S. Xie, PDRF: A general dispersion relation solver for 
%        magnetized multi-fluid plasma, Computer Physics Communications,
%        185, 670-675 (2014). http://dx.doi.org/10.1016/j.cpc.2013.10.012
%
% test run:
%   close all; clear; clc; 
%   [M,A,w1,w2,kk,ww1,ww2,betasz,vA,csz2,c2,rhos0,wcs,qs,ms,wps,B0]=pdrf_SI;
%
% Note: 
%   Although this solver has been verified under some basic benchmarks,
%   you'd better verify it before use it for more complicated problems.
%   See [Xie2014] for details.
% 
% 2016-04-14 18:56 ver 1.1, fix bugs for:
%        1. roundoff error of eig(), using vpa() sym format
%        2. P\neq0 term, corresponding matrix elements rho_s0 -> n_s0
% Ackn.: Nicholas Francken and Yana Maneva at the KU Leuven in Belgium
%        for the benchmarks report, which is helpful to fix these two bugs.
%
% URL: http://hsxie.me/codes/pdrf/
%

% run pdrk;
% global tmpcwp tmpwcs1; tmpcwp=cwp; tmpwcs1=wcs1;
function [M,A,w1,w2,kk,ww1,ww2,betasz,vA,csz2,c2,rhos0,...
                              wcs,qs,ms,wps,B0]=pdrf_SI  % some test runs

%     close all; clear; clc;

    global s qs ms ns0 vs0x vs0y vs0z csz2 csp2 epsnx epsny rhos0 Q0 ...
        J0x J0y J0z nu B0 c2 mu0 epsilon0 gammaTs Psz Psp ...
        betasz betasp Deltas vA wcs wps;
    global tmpcwp tmpwcs1;
    
    initial; % call initial

    rel=0;  % rel=1, relativistic; else non-relativistic.

    method=2; % 1: single k; 2: w(k); 3: w(theta); 4: w(kx,kz);

    % Here are some default test runs with output data or figures.
    % Modify here for your new cases.
    if (method>=3)
        h=figure('unit','normalized','Position',[0.01 0.27 0.6 0.5]);
        set(gcf,'DefaultAxesFontSize',15);
    end
    
    if(method==2)        
%         kk=[0.0:1e-7:1e-5,2e-5:5e-6:2e-4,3e-4:1e-4:6e-3];
        kk=0.0:1e-5:1e-3;
        theta=60*pi/180;
        nk=length(kk);
        ww=zeros(nk,4*s+6); ww1=ww; ww2=ww;
        for ik=1:nk
            kx=kk(ik)*sin(theta);
            kz=kk(ik)*cos(theta);
            [M,A,w1,w2]=pdrfsolver(kx,kz,rel);
            ww1(ik,:)=w1;
            ww2(ik,:)=w2;
        end
        ns=4*s+6;
        for is=1:ns
            clr=colormap(lines);
            clr=clr(floor(is/ns*length(clr)),:);
            subplot(121);
            plot(kk*tmpcwp,real(ww1(:,is))/tmpwcs1,'Color',clr,'LineWidth',2); 
            hold on; 
%             xlabel('k'); ylabel('\omega_r'); xlim([min(kk),max(kk)]);
%             title(['(a) \omega_r vs k, \theta=',...
%                 num2str(theta*180/pi),'^{\circ}']);
            subplot(122); 
            plot(kk*tmpcwp,imag(ww2(:,is))/tmpwcs1,'.','Color',clr,'LineWidth',2); 
            hold on; 
%             xlabel('k'); ylabel('\omega_i'); xlim([min(kk),max(kk)]);
%             title(['(b) \omega_i vs k, m_i/m_e=',num2str(ms(2)/ms(1))]);
        end
        betasz
        betasp
    elseif(method==3)
        
        k=2.0; dthe=pi/100;
        tt=0:dthe:pi/2;
        nt=length(tt);
        ww=zeros(nt,4*s+6); ww1=ww; ww2=ww;
        for it=1:nt
            kx=k*sin(tt(it));
            kz=k*cos(tt(it));
            [M,A,w1,w2]=pdrfsolver(kx,kz,rel);
            ww1(it,:)=w1;
            ww2(it,:)=w2;
        end
        ns=2*s+1;
        for is=1:ns
            clr=colormap(lines);
            clr=clr(floor(is/ns*length(clr)),:);
            subplot(121);
            plot(tt,real(ww1(:,is)),'x','Color',clr,'LineWidth',2); 
            hold on; xlabel('\theta'); ylabel('\omega_r');
            title(['(c) \omega_r vs \theta, k=',num2str(k)]); 
            xlim([min(tt),max(tt)]);
            subplot(122); 
            plot(tt,imag(ww2(:,is)),'.','Color',clr,'LineWidth',2); 
            hold on; xlabel('\theta'); ylabel('\omega_i'); 
            xlim([min(tt),max(tt)]);
            title('(d) \omega_i vs \theta');
        end
        
    elseif(method==4)
        
        [kkx,kkz]=meshgrid(0.01:0.04:3.0,0.01:0.04:2.0);
        [row,col]=size(kkx);
        ww=zeros(row,col,4*s+6);
        ww1=ww;
        ww2=ww;

        for ir=1:row
            for ic=1:col
                kx=kkx(ir,ic);
                kz=kkz(ir,ic);
                [M,A,w1,w2]=pdrfsolver(kx,kz,rel);
                ww1(ir,ic,:)=w1;
                ww2(ir,ic,:)=w2;
            end
        end
        ns=4*s+6;
        subplot(121);
        for is=1:ns
            surf(kkx,kkz,real(ww1(:,:,is))); hold on;
        end
        xlabel('k_x'); ylabel('k_z'); 
        zlabel('\omega_r');

        subplot(122);
        for is=1:ns
            surf(kkx,kkz,imag(ww2(:,:,is))); hold on;
        end
        xlabel('k_x'); ylabel('k_z'); 
        zlabel('\omega_i');  
        
    else % single k, note: k can also be complex number
        
        k=0.1; theta=pi/3; 
        kx=k*sin(theta); kz=k*cos(theta);
%         kx=0.2; kz=0.01;
        [M,A,w1,w2]=pdrfsolver(kx,kz,rel);
        w=w1  % output
        
    end
    
%     print('-dpng','pdrf.png');

end

function initial

    global s qs ms ns0 vs0x vs0y vs0z csz2 csp2 epsnx epsny rhos0 Q0 ...
        J0x J0y J0z nu B0 c2 mu0 epsilon0 gammaTs Psz Psp ...
        betasz betasp Deltas vA;
    
    % read input parameters
    par = importdata('pdrf_SI.in', ' ', 1); 
    [s, col]=size(par.data);
    if(col~=10)
        disp('Wrong input data !!!');
    end

    for j=1:s
        qs(j)=par.data(j,1);
        ms(j)=par.data(j,2);
        ns0(j)=par.data(j,3);
        vs0x(j)=par.data(j,4);
        vs0y(j)=par.data(j,5);
        vs0z(j)=par.data(j,6);
        csz2(j)=par.data(j,7)^2;
        csp2(j)=par.data(j,8)^2;
        epsnx(j)=par.data(j,9);
        epsny(j)=par.data(j,10);
    end
    rhos0=ns0.*ms;
    Q0=sum(qs.*ns0);
    J0x=sum(qs.*ns0.*vs0x);
    J0y=sum(qs.*ns0.*vs0y);
    J0z=sum(qs.*ns0.*vs0z);
    if(Q0~=0)
        disp('Non-neutral charged !!!');
    end

    % collision parameters, modify them as you like
    nue=0.0; nuee=0.0; nui=0.0; nuii=0; nuei=0.0;
    nu=zeros(s,s);
    for i=1:s
        if(ms(i)<=1)
            nu(i,i)=nue;
        else
            nu(i,i)=nui;
        end
        for j=(i+1):s
            if(ms(i)<=1 && ms(j)<=1)
                nu(i,j)=nuee;
                nu(j,i)=nuee;
            elseif(ms(i)<=1 && ms(j)>1)
                nu(i,j)=nuei;
                nu(j,i)=nuei;
            else
                nu(i,j)=nuii;
                nu(j,i)=nuii;
            end
        end
    end

    % other initial parameters, SI unit
    B0=100.0E-9; c2=(2.9979E8)^2; epsilon0=8.854E-12 ; mu0=1/(c2*epsilon0);
    
    gammaTs=0.*qs+5/3; % adiabatic coefficients, modify here for para/perp
    Psz=csz2.*rhos0./gammaTs;
    Psp=csp2.*rhos0./gammaTs;
    if (B0~=0)
        betasz=2*mu0.*Psz./B0^2;
        betasp=2*mu0.*Psp./B0^2;
        vA=B0./sqrt(mu0.*rhos0); % 

        Deltas=-(Psp-Psz)/B0;
    else
        Deltas=0.*qs;
    end
    
end

function [M,A,w1,w2,X1,X2]=pdrfsolver(kx,kz,rel)
% main solver of pdrf, call as
%  (1). [M,A,w1,w2,X1,X2]=pdrfsolver(kx,kz,rel)
%          for both eigenvalues and polarizations
%  (2). [M,A,w1,w2]=pdrfsolver(kx,kz,rel)
%          for just eigenvalues
% (1) requires more run time than (2)

    global s qs ms ns0 vs0x vs0y vs0z csz2 csp2 epsnx epsny rhos0 Q0 ...
        J0x J0y J0z nu B0 c2 mu0 epsilon0 gammaTs Psz Psp ...
        betasz betasp Deltas vA wcs wps;

    % eigen matrix
    M=zeros(4*s+6,4*s+6);
    A=eye(4*s+6,4*s+6);

    wcs=qs.*B0./ms;
    wps=sqrt(ns0.*qs.^2./epsilon0./ms);
    ind2=4*s;
    for j=1:s
        ind=(j-1)*4;

        % dn~n
        M(1+ind,1+ind)=-1i*(kx*vs0x(j)+kz*vs0z(j));
        % dn~v
        M(1+ind,2+ind)=-1i*kx*ns0(j)-epsnx(j)*ns0(j);
        M(1+ind,3+ind)=-epsny(j)*ns0(j);
        M(1+ind,4+ind)=-1i*kz*ns0(j);

        ajpq=eye(3); % ajpq for A
        if(rel==1)
            vj0=sqrt(vs0x(j)^2+vs0y(j)^2+vs0z(j)^2);
            gammaj0=1/sqrt(1-vj0^2/c2);
            ajpq(1,1)=gammaj0+gammaj0^3*vs0x(j)^2/c2;
            ajpq(1,2)=gammaj0^3*vs0x(j)*vs0y(j)/c2;
            ajpq(1,3)=gammaj0^3*vs0x(j)*vs0z(j)/c2;
            ajpq(2,1)=gammaj0^3*vs0x(j)*vs0y(j)/c2;
            ajpq(2,2)=gammaj0+gammaj0^3*vs0y(j)^2/c2;
            ajpq(2,3)=gammaj0^3*vs0y(j)*vs0z(j)/c2;
            ajpq(3,1)=gammaj0^3*vs0x(j)*vs0z(j)/c2;
            ajpq(3,2)=gammaj0^3*vs0z(j)*vs0y(j)/c2;
            ajpq(3,3)=gammaj0+gammaj0^3*vs0z(j)^2/c2;
        end
        A((2:4)+ind,(2:4)+ind)=ajpq;

        bjpq=-1i*(kx*vs0x(j)+kz*vs0z(j)).*ajpq; % bjpq for M, nu(j,j) below
        M((2:4)+ind,(2:4)+ind)=bjpq;

        % dv ~ n & v
%         M(2+ind,1+ind)=-1i*kx*csp2(j)/rhos0(j); % wrong
%         M(4+ind,1+ind)=-1i*kz*csz2(j)/rhos0(j);
        M(2+ind,1+ind)=-1i*kx*csp2(j)/ns0(j); % 16-04-13 12:34, bug fixed
        M(4+ind,1+ind)=-1i*kz*csz2(j)/ns0(j);
        
        M(2+ind,3+ind)=M(2+ind,3+ind)+wcs(j);
        M(3+ind,2+ind)=M(3+ind,2+ind)-wcs(j);

        % dv~E
        M(2+ind,1+ind2)=qs(j)/ms(j);
        M(3+ind,2+ind2)=qs(j)/ms(j);
        M(4+ind,3+ind2)=qs(j)/ms(j);    
        % dv~B
        M(2+ind,5+ind2)=-qs(j)/ms(j)*vs0z(j);
        M(2+ind,6+ind2)=qs(j)/ms(j)*vs0y(j);
        M(3+ind,4+ind2)=qs(j)/ms(j)*vs0z(j);
        M(3+ind,6+ind2)=-qs(j)/ms(j)*vs0x(j);
        M(4+ind,4+ind2)=-qs(j)/ms(j)*vs0y(j);
        M(4+ind,5+ind2)=qs(j)/ms(j)*vs0x(j);
        
        M(2+ind,4+ind2)=-1i*kz*Deltas(j)/rhos0(j);
        M(3+ind,5+ind2)=-1i*kz*Deltas(j)/rhos0(j);
%         M(4+ind,4+ind2)=M(4+ind,4+ind2)-...
%             1i*kx*1*Deltas(j)/rhos0(j); % without factor "2" correction
        M(4+ind,4+ind2)=M(4+ind,4+ind2)-...
            1i*kx*2*Deltas(j)/rhos0(j); % with "2" correction, mirror mode

        % dE~n
        M(1+ind2,1+ind)=-qs(j)*vs0x(j)/epsilon0;
        M(2+ind2,1+ind)=-qs(j)*vs0y(j)/epsilon0;
        M(3+ind2,1+ind)=-qs(j)*vs0z(j)/epsilon0;
        
        % dE~v
        M(1+ind2,2+ind)=-qs(j)*ns0(j)/epsilon0;
        M(2+ind2,3+ind)=-qs(j)*ns0(j)/epsilon0;
        M(3+ind2,4+ind)=-qs(j)*ns0(j)/epsilon0;

    end

    % dE~B
    M(1+ind2,5+ind2)=-1i*kz*c2;
    M(2+ind2,4+ind2)=1i*kz*c2;
    M(2+ind2,6+ind2)=-1i*kx*c2;
    M(3+ind2,5+ind2)=1i*kx*c2;

    % dB~E
    M(4+ind2,2+ind2)=1i*kz;
    M(5+ind2,1+ind2)=-1i*kz;
    M(5+ind2,3+ind2)=1i*kx;
    M(6+ind2,2+ind2)=-1i*kx;

    % set collision nu(i,j) terms to M
    for i=1:s
        for j=1:s
            col=(i-1)*4;
            row=(j-1)*4;
            M(2+row,2+col)=M(2+row,2+col)+nu(i,j);
            M(3+row,3+col)=M(3+row,3+col)+nu(i,j);
            M(4+row,4+col)=M(4+row,4+col)+nu(i,j);
        end
    end

    % solution
    oldmethod=2; % default oldmethod=2, less roundoff error
    % oldmethod=1; faster, but may less accurate
    if(nargout==6) % give also eigenvector, for polarizations        
        if(oldmethod==1)
            [V,D]=eig(M,A);
            d=diag(D);
        else % to keep all eigenvalue be accurate to digits(n)
            MA=A\M;
            [V0,D0]=eig(vpa(MA,16));  % sym format calculation, less roundoff error
            V=double(V0); d0=diag(D0);
            d=double(d0); % convert sym format to double
        end
    else
        if(oldmethod==1)
            d=eig(M,A); % eigenvalue, i.e., dispersion relation solutions
        else % better for max(M_ij)/min(M_ij)>10^16,  16-04-12 20:41
            % to keep all eigenvalue be accurate to digits(n)
            MA=A\M; % A^{-1}*M, i.e., M*X=lambda*A*X -> A^{-1}*M*X=lambda*X
            d0=vpa(eig(MA),16);  % sym format calculation, less roundoff error
%             d0=eig(vpa(MA,16)); % seems slow
            d=double(d0); % convert sym format to double
        end
    end
    
    wtmp=1i*d;
    [wr,inw1]=sort(real(wtmp),'descend');
    [wi,inw2]=sort(imag(wtmp),'descend');
    w1=wtmp(inw1); % sort by real part
    w2=wtmp(inw2); % sort by imag part
    
    if(nargout==6)
        X1=V(:,inw1); % sort by real part
        X2=V(:,inw2); % sort by imag part
    end
    
end
