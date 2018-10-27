% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2014-06-01 17:00
% pdrk_em3d.m, Plasma Dispersion Relastion solver, kinetic, EM3D,
% bi-Maxwellian equilibrium distribution with parallel drift.
% Transform to matrix eigenvalue problem lambda*X=M*X.
% J-pole approximation for Z(zeta)=sum(b_j/(zeta-c_j))
%
% sparse, eigs
%
% 2014-08-22 01:09, support T_perp/T_para \neq 1,  vds \neq 0
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
% 18-10-13 10:37 update to with loss-cone distribution, which thus can
% handle all the same input cases as WHAMP [Ronnmark1982] code
% ------------------------------------------------------------------------
% This file is the kernel part of pdrk-em3d code.
% Modify sp and nw if you only need one or several solutions around the
% initial guess value wg.

wwp=zeros(npa,npb,nw);
[ppa,ppb]=ndgrid(pa,pb);
kxx=0.*ppa; % for store kx
kzz=0.*ppa; % for store kz
kk=0.*ppa; % for store k
tt=0.*ppa; % for store theta
%tic;
for jpa = 1:npa % scan 1st parameter a
  
  for jpb = 1:npb % scan 2nd parameter b
    
    % ------- to update 2018-10-18 20:10 --------
    if(iloga==0)
      par(ipa)=pa(jpa);
    else
      par(ipa)=10^(pa(jpa));
    end
    if(ipa~=ipb)
      if(ilogb==0)
        par(ipb)=pb(jpb);
      else
        par(ipb)=10^(pb(jpb));  
      end
    end
    
    k=abs(par(1))/cwp; % k only >=0
    theta=par(2); % the angle between k and B0
    if(ipa>2 && ipb>2) % !! to update
      kz=par(3)/cwp;
      kx=par(4)/cwp;
    else
      kz=cos(theta*pi/180)*k;
      kx=sin(theta*pi/180)*k;
    end
    k=sqrt(kx*kx+kz*kz); % update k
    theta=angle(kz+1i*kx)*180/pi; % update theta
    % !! if you want scan other parameters, modify here, e.g.,
    % betasz=par(5); % or vds(1)=par(5); % etc
    % ---------------
    
    if(iem==1) % electromagnetic run
      run pdrk_em3d_matrix;
    else % electrostatic run, 18-10-20 20:02
      run pdrk_es3d_matrix;
    end
    
    % calculate the solution use either eig() or eigs()
    if(sp==0) % solve all the roots, slow for large N
      % d=eig(full(M));
      % it's found in pdrf (fluid version) that vpa(eig) may more accurate
      d0=vpa(eig(full(M)),16);d=double(d0);
    else % solve only nw solutions, fast for large N, similar to WHAMP
      % d=eigs(M,nw,wg);
      if(iem==1 && iout==2 && icalp==1) % for pdrk_output.m
        disp(['jpa=',num2str(jpa),', jpb=',num2str(jpb)]);
        eps=1i*1e-8; % 18-10-06 10:29
        wg=wws(jpa,jpb,jpl)*wcs1+eps; % add eps to avoid singular warning
        %d0=vpa(eigs(M,nw,wg),16);d=double(d0);
        [V,D]=eigs(M,nw,wg);
        %[V0,D0]=vpa(eigs(M,nw,wg)); V=double(V0); D=double(D0);
        dEB=V((NN-5):NN,1); % polarization data (dE,dB)
        dEB=dEB/dEB(1); % 18-10-19 17:04
%         ctmp=sqrt(real(dEB(1))^2+real(dEB(2))^2+real(dEB(3))^2)+...
%             1i*sqrt(imag(dEB(1))^2+imag(dEB(2))^2+imag(dEB(3))^2);
        ctmp=sqrt(real(dEB(1))^2+real(dEB(2))^2+real(dEB(3))^2)+...
            imag(dEB(1))^2+imag(dEB(2))^2+imag(dEB(3))^2;
        dEB=dEB/ctmp; % normalized (E,B) % 18-10-18 07:14
        Pola(jpa,jpb,jpl,1:6)=dEB; % store dEB
        wws2(jpa,jpb,jpl)=D(1,1)/wcs1; % store new wws, in case
        
        dEx=dEB(1);
        dEy=dEB(2);
        dEz=dEB(3);
        dBx=dEB(4);
        dBy=dEB(5);
        dBz=dEB(6);
        UE=(dEx*conj(dEx)+dEy*conj(dEy)+dEz*conj(dEz))*epsilon0*0.5; % electric energy
        UB=(dBx*conj(dBx)+dBy*conj(dBy)+dBz*conj(dBz))/mu0*0.5; % magnetic energy
        Pola(jpa,jpb,jpl,7)=UE;
        Pola(jpa,jpb,jpl,8)=UB;
        
        % kdotE=abs(kx*dEx+kz*dEz); % to add more
        % kcrossE=abs(sqrt((-kz*dEy)^2+(kz*dEx-kx*dEz)^2+(kx*dEy)^2)); %
        % abs(k \dot E), abs(k \times E)
        % Pola(jpa,jpb,jpl,9)=kdotE; %
        % Pola(jpa,jpb,jpl,10)=kcrossE; %
        
      else % calcute also the polarizations
        d0=vpa(eigs(M,nw,wg),16);d=double(d0);
      end
    end
    omega=d;
    [wi,ind]=sort(imag(omega),'descend'); % sort the solution by growth rate
    w=omega(ind);
    wg=w(1); % update the initial guess for pa(jpa,jpb+1) use previous solution
    
    if(jpb==1) % initial guess for pb(jpa+1,jpb)
      wg0=wg;
    end

    wwp(jpa,jpb,1:nw)=w(1:nw);
    
    if(icalp==0) % store k, theta in the first run
      kxx(jpa,jpb)=kx;
      kzz(jpa,jpb)=kz;
      kk(jpa,jpb)=k;
      tt(jpa,jpb)=theta;
    end
  end %jpb
  
  wg=wg0; % update the new initial guess
  
end %jpa
%runtime=toc;
if(icalp==0)
  ww=wwp; % store all the solutions in the first run
end
