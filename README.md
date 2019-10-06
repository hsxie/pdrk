PDRK (http://hsxie.me/codes/pdrk/, https://github.com/hsxie/pdrk) and PDRF (http://hsxie.me/codes/pdrf/, https://github.com/hsxie/pdrf) have been renamed and updated to a new code BO (‘波’, i.e., 'wave' in Chinese), which is a state-of-art plasma wave and instability analysis tool and is being used by more and more researchers. The new version is not included in this github page yet and can be found at http://code.ennresearch.com/bo/, or http://dx.doi.org/10.17632/cvbrftzfy5.1.

Ref: H. S. Xie, BO: A unified tool for plasma waves and instabilities analysis, Comp. Phys. Comm., 244, 343 (2019), https://doi.org/10.1016/j.cpc.2019.06.014.

![Image text](https://github.com/hsxie/pdrk/blob/master/bo_poster_191006.jpg)

Welcome to share this poster: https://github.com/hsxie/pdrk/blob/master/bo_poster.pdf!

Hua-sheng XIE, huashengxie[AT]gmail.com, xiehuasheng[AT]enn.cn

CCF-ENN, Langfang, China

2019-10-06 08:45

-------------------------------------------------------------------------

PDRK is a general and powerful Kinetic Dispersion Relation solver for 
magnetized Plasma, which is the first kinetic dispersion relation solver 
that can give all the important solutions at one time without requiring
initial guess for root finding.

This is an update version of pdrk-em3d code, which have fixed the bugs in 
the original version in [Xie2016] paper.

Also, this version is more user friendly. Several sub-files are added to 
separate the structures of the code.

And, most importantly, this new version we have overcome the major drawback
in the original version of pdrk.

% This version (2018-10-19 version): 

%  1. With loss-cone, drift bi-Maxwellian distribution;

%  2. With module to separate dispersion surfaces;

%  3. Support output polarizations, etc;

%  4. Support EM3D, ES3D and ES1D.

-------------------------------------------------------------------------

.Readme.txt            % this file

./doc

 - *.pdf               % documents, equations, user guide


./examples             % some typical test examples

 - * 


./code

 - pdrk_main.m         % use this file to run the code


 ./input
 
  - pdrk.in             % input parameters for each species
  
  - pdrk_setup.m        % initial setup of how to scan the parameters
  
  - pdrk_wpdat.m        % initial data of omega for select plots


 ./modules
 
  - pdrk_initialize.m   % initialize the parameters for run the code
  
  - pdrk_kernel.m       % kernel, most important part of pdrk code
  
  - pdrk_em3d_matrix.m   % set the pdrk em3d matrix, used by kernel
  
  - pdrk_es3d_matrix.m   % set the pdrk es3d matrix, used by kernel
  
  - pdrk_plot_all.m     % output data and plot
  
  - pdrk_plot_select.m  % subtle file for separate dispersion surfaces
  
  - pdrk_output.m       % calculate the polarizations and output data
  

 ./output               % directory for store figures and data
 
  - * 


-------------------------------------------------------------------------

I will long time support this code. So, please feedback to me of your
suggestions/experience.

If you use this code, please cite:

[Xie2016] Huasheng Xie and Yong Xiao, PDRK: A General Kinetic Dispersion 
Relation Solver for Magnetized Plasma, Plasma Science and Technology, 18, 
2, 97 (2016). Update/Bugs fixed at http://hsxie.me/codes/pdrk/ or 
https://github.com/hsxie/pdrk.

Enjoy!

Hua-sheng XIE, huashengxie@gmail.com

FRI-ENN, China

2018-10-27 17:31


-------------------------------------------------------------------------

------------------------------ PDRK history -----------------------------

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


% 18-10-03 13:14 This version should be bugs free now.

% 18-10-13 10:37 update to with loss-cone distribution, which thus can

% handle all the same input cases as WHAMP [Ronnmark1982] code.

% 18-10-20 18:14 update to support electrostatic loss-cone.

