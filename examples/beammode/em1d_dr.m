% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2012-01-10 14:01
% On Numerical Calculation of the Plasma Dispersion Function
% Ion acoustic wave, ES1D
% 2018-10-15 15:46 Modified to EM1D (only kpara)
% clear;clc;

zeta=@(x)faddeeva(x)*1i*sqrt(pi);
% f=@(x,k)k*k+(1+x*zeta(x))+1/Ti*(1+(x/sqrt(Ti/mi))*zeta(x/sqrt(Ti/mi)));

pm=1; % pm=1 and -1 for different branches
f=@(x,k)x^2-k*k*c2...
    +wps2(1)*(x-k*vds(1))/(k*vtzs(1))*zeta((x-k*vds(1)+pm*wcs(1))/(k*vtzs(1)))...
    +wps2(2)*(x-k*vds(2))/(k*vtzs(2))*zeta((x-k*vds(2)+pm*wcs(2))/(k*vtzs(2)))...
    +wps2(3)*(x-k*vds(3))/(k*vtzs(3))*zeta((x-k*vds(3)+pm*wcs(3))/(k*vtzs(3)))...
    +wps2(4)*(x-k*vds(4))/(k*vtzs(4))*zeta((x-k*vds(4)+pm*wcs(4))/(k*vtzs(4)));
w=[];

kmin=0.01;dk=0.01;kmax=0.5;
% k=(kmax:-dk:kmin)/cwp;
k=(kmin:dk:kmax)/cwp;
% x0=(1.0-0.0i)*sqrt(2)*kmin;
% x0=1.2-0.5i;
x0=(0.05+0.04i)*wcs1;
% x0=(0.005+0.0i)*wcs1;
for kk=k
    options=optimset('Display','off');
    x=fsolve(f,x0,options,kk);
    x0=x;
    w=[w,x];
end
wre=real(w);wie=imag(w);
%%
% figure;set(gcf,'DefaultAxesFontSize',15);
% plot(k*cwp,wre/wci,'r+',k*cwp,wie/wci,'o','LineWidth',2);
subplot(121); hold on;
plot(k*cwp,wre/wcs1,'k+','LineWidth',2);
subplot(122); hold on;
plot(k*cwp,wie/wcs1,'k+','LineWidth',2);
grid on;
% xlabel('kc/\omega_{pi}');ylabel('\omega/\omega_{cp}');
% xlim([kmin,kmax]); title('EM1D \theta=0');
%%
if(pm==-1)
subplot(121);
legend('pdrk-1','pdrk-2','pdrk-3','pdrk-4','EM1D DR-1','EM1D DR-2',...
    'location','northwest'); legend('boxoff');
end
% legend('\omega_r','\gamma','location','northwest'); legend('boxoff');
