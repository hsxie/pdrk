% 18-10-05 08:00 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO),
% etc ...
% 18-10-13 10:37 update to with loss-cone distribution, which thus can
% handle all the same input cases as WHAMP [Ronnmark1982] code.
% This version: 
%   1. With loss-cone, drift bi-Maxwellian distribution;
%   2. With module to separate dispersion surfaces;
%   3. Support output polarizations, etc;
%   4. Support EM3D, ES3D and ES1D.
% =======================================================================
% Use this file to start run pdrk.

close all;
clear; clc;

% % set initial parameters in 'pdrk.in' and 'pdrk_setup.m'
disp('------------------- PDRK version 2018-10-27 -------------------');
tic; runtime=0;
disp('run ./input/pdrk_setup.m ...');
run ./input/pdrk_setup;
runtime1=toc; runtime=runtime+runtime1;
disp([' use ',num2str(runtime1),' s, toal ',num2str(runtime),'s']);

tic;
% % initialize the parameters for 'pdrk_kernel.m'
disp('run ./modules/pdrk_initialize.m ...');
run ./modules/pdrk_initialize;
runtime2=toc; runtime=runtime+runtime2;
disp([' use ',num2str(runtime2),' s, toal ',num2str(runtime),'s']);
if(ireturn==1)  
  disp(['(ipa,ipb) = (',num2str(ipa),', ',num2str(ipb),...
    '), [see ''setup.m'', 1: k; 2: theta; 3: kz; 4: kx]']);
  disp(strscan);
  return;
end

% % print some basic info
disp('---------------------------------------------------------------');
if(iem==1)
  disp(['iem = ',num2str(iem),', this is an * electromagnetic * run.']);
else
  disp(['iem = ',num2str(iem),', this is an * electrostatic * run.']);
end
disp(['N = ',num2str(N),', [Number harmonics. ! Try larger N to',...
    ' make sure the results are convergent]']);
disp(['J = ',num2str(J),', [J-pole, usually J=8 is sufficient; ',...
    'other choice: 4 (faster), 12 (more accurate)]']);
disp(['sp = ',num2str(sp),', [sp=1, sparse eigs() to obtain nw0 ',...
    'solutions; sp=0, eig() to obtain all solutions]']);

disp(['nw0 = ',num2str(nw0),', [solve nw0 solutions one time]']);
disp(['iout = ',num2str(iout),', [=1, only (omega,gamma); =2, ',...
    'also calculate (E,B)]']);

disp(['(npa,npb) = (',num2str(npa),', ',num2str(npb),...
    '), [number scan for the (1st, 2nd) variables]']);
disp(['(ipa,ipb) = (',num2str(ipa),', ',num2str(ipb),...
    '), [see ''setup.m'', 1: k; 2: theta; 3: kz; 4: kx]']);
disp(['(iloga,ilogb) = (',num2str(iloga),', ',num2str(ilogb),...
    '), [=0, linear scan; =1, 10^(pa or pb) scan]']);
disp(['pa=',num2str(min(pa)),':',num2str((pa(end)-pa(1))/(npa-1+1e-10)),...
    ':',num2str(max(pa)),', [1st parameter range]']);
disp(['pb=',num2str(min(pb)),':',num2str((pb(end)-pb(1))/(npb-1+1e-10)),...
    ':',num2str(max(pb)),', [2nd parameter range]']);
disp(['This run: pa = ',strpa,', pb = ',strpb,', ',strscan]);

disp('---------------------------------------------------------------');
disp(['S [number of species] = ',num2str(S)]);
disp(['B0 [backgroud magnetic field, Tesla] = ',num2str(B0)]);

disp('---------------------------------------------------------------');

disp(['qs0 [charge, q/e] = ',num2str(qs0)]);
disp(['ms0 [mass, m/mp] = ',num2str(ms0)]);
disp(['ns0 [desity, m^-3] = ',num2str(ns0)]);
disp(['Tzs0 [parallel temperature, eV] = ',num2str(Tzs0)]);
disp(['Tps0 [perp temperature, eV] = ',num2str(Tps0)]);
disp(['alphas [loss-cone size/anisotropic] = ',num2str(alphas)]);
disp(['Deltas [loss-cone depth, =0 (max) to 1 (no)] = ',num2str(Deltas)]);
disp(['vds0 [para drift velocity, vds/c] = ',num2str(vds0)]);
disp(['rsa [core ratio] = ',num2str(rsab(1,:))]);
disp(['rsb [loss cone ratio] = ',num2str(rsab(2,:))]);

disp('---------------------------------------------------------------');

disp(['lambdaDs [Debye length, m] = ',num2str(lambdaDs)]);
disp(['wps [plasma frequency, Hz] = ',num2str(wps)]);
disp(['wcs [cyclotron frequency, Hz] = ',num2str(wcs)]);

disp(['rhocs [cyclotron radius, m] = ',num2str(rhocs)]);
disp(['wps [plasma frequency, Hz] = ',num2str(wps)]);
disp(['lmdT [Tpara/Tperp] = ',num2str(lmdT)]);
disp(['betasz [parallel beta] = ',num2str(betasz)]);
disp(['betasp [perp beta] = ',num2str(betasp)]);
disp(['vA [=B0/sqrt(mu0*sum(ms.*ns0)), Alfven speed, m/s] = ',num2str(vA)]);
disp(['c [speed of light, m/s] = ',num2str(sqrt(c2))]);

disp('====- In PDRK plot/output, k -> k*cwp, omega -> omega/wcs1 -====');
disp(['---- ! Set the 1st species to be ion in ''pdrk.in'', if',10,...
    '---- ! you hope wcs1=omega_ci and cwp=c/omega_pi.']);

disp(['wcs1 [1st species cyclotron frequency, Hz] = ',num2str(wcs1)]);
disp(['wps1 [1st species plasma frequency, Hz] = ',num2str(wps1)]);
disp(['cwp [c/wps1, m] = ',num2str(cwp)]);

disp('---------------------------------------------------------------');

%%
tic;
% the kernel part of pdrk_em3d code, donot need change for most cases
% Or, only modify 'pdrk_kernel.m' lines [43-65, xxx] %
icalp=0; % do not calculate polarization in the first step
disp('run ./modules/pdrk_kernel.m ...');
run ./modules/pdrk_kernel;
runtime3=toc; runtime=runtime+runtime3;
disp([' use ',num2str(runtime3),' s, toal ',num2str(runtime),'s']);

tic;
% plot all the solutions
disp('run ./modules/pdrk_plot_all.m ...');
run ./modules/pdrk_plot_all;
runtime4=toc; runtime=runtime+runtime4;
disp([' use ',num2str(runtime4),' s, toal ',num2str(runtime),'s']);
disp('---------------------------------------------------------------');

% % plot select dispersion surface, very subtle and can be further improve
input(['Please update ''wpdat'' in ''./input/pdrk_wpdat.m'' firstly. ',10,...
    'This step is to determine which branch(es) you hope to output/store.',...
    10,'After set the ''wpdat'', press any key to continue. ']);
disp('---------------------------------------------------------------');
tic;
disp('run ./modules/pdrk_plot_select.m use ./input/pdrk_wpdat ...');
% run ./modules/pdrk_plot_select;
run ./input/pdrk_wpdat; % 18-10-21 18:04
runtime5=toc; runtime=runtime+runtime5;
disp([' use ',num2str(runtime5),' s, toal ',num2str(runtime),'s']);

input('Press any key to continue for the last step: run output.');
disp('---------------------------------------------------------------');
tic;
% % output the omega and/or polarization results to data file
disp('run ./modules/pdrk_output.m ...');
run ./modules/pdrk_output;
runtime6=toc; runtime=runtime+runtime6;
disp([' use ',num2str(runtime6),' s, toal ',num2str(runtime),'s']);
disp(['Finished! Figures/data have been saved to ',savepath,'.']);
