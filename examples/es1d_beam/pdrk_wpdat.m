% 18-10-19 17:56 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO),
% etc ...
% Initial data for run pdrk_plot_select

% Search the most close dispersion surfaces to these data.
% Initial data for find the corresponding dispersion surfaces.
% Please use pdrk_plot_all.m to visualize all the solutions, and then
% modify here the initial point of which mode(s) you want plot/store.

% wpdat(:,1) is pa; wpdat(:,2) is pb for 2D scan and arbitrary for 1D scan;
% wpdat(:,3) is Re or Im(omega)

% rex=abs(lambdaDs(1)/cwp); % rescale the x-axis, 18-10-21 08:54
rey=abs(lambdaDs(1)/cwp); % rescale the y-axis
rez=abs(wcs(1)/sqrt(sum(wps2))); % rescale the omega-axis, 18-10-27 10:05
rex=abs(sqrt(1/sum(1./lambdaDs.^2))/cwp);

wpdat=[151,0,0.4451i;
    171,0,-0.2437i;
    231,0,-0.04389i;
    150,0,-0.5583i;
  ];

% wpdat=[
%    311,0,1.425;
%    511,0,2.31;
%    511,0,3.282;
%    1171,0,4.172;
%   ];

% wpdat=[361,0,0.0011i;
%    201,0,0.0009i;
%    301,0,3.58e-4*1i;
%    511,0,-3.79e-4*1i;
%   ];
  
% wpdat=[311,0,-0.209i;
%   ];
  
% wpdat=[65,0,5.28e-4*1i;
%     73,0,7.5e-4*1i;
%   32,0,2.447;
%   9,0,3.46;
%   24,0,4.236;
%   ];


run ../modules/pdrk_plot_select;
%%
% calculate omega of k*lambdaD=0.2
kinterp=0.2;
winterp=interp1(rex*pas,rez*wws(:,1,1),kinterp);
subplot(122);
xlim([0,0.5]);
% ylim([-1e-3,2e-3]);
plot(kinterp,imag(winterp),'x'); hold on;
text(kinterp+0.02,imag(winterp),[num2str(kinterp),',',num2str(imag(winterp),4)]);

subplot(121);
xlim([0,0.5]);
plot(kinterp,real(winterp),'x'); hold on;
text(kinterp+0.02,real(winterp),[num2str(kinterp),',',num2str(real(winterp),4)]);

