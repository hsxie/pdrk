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
% rey=abs(lambdaDs(1)/cwp); % rescale the y-axis
% rez=abs(wcs(1)/sqrt(sum(wps2))); % rescale the omega-axis, 18-10-27 10:05
rex=1;
rey=1;
rez=1;


wpdat=[%70,0,-152i; % for 2D scan
  90.5,85,-261.5i;
  %94.5,85,-433i;
  ];

%%
run ../modules/pdrk_plot_select;
% subplot(122);xlim([0,2]);ylim([-1e-3,2e-3]);subplot(121);xlim([0,2]);
