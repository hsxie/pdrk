% 18-10-05 08:00 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO),
% etc ...
% Setup for run pdrk-em3d
% 18-10-13 09:52 update to with loss-cone distribution

% savepath='../output/'; % choose where to save the figures and outputs
savepath='../examples/es1d_beam/';
iem=0; % =1, electromagnetic run; =0, electrostatic run

N=1; % Number harmonics, enlarge N to make sure the results are convergent
J=8; % J-pole, usually J=8 is sufficient; other choice: 4, 12
sp=0; % sparse eigs(), sp=1; or eig() for all solution, sp=0
wg0=1.5+0.01i; % inital guess for sp=1, normalized by omega_{c1}; not required if sp=0

B0=143.5E-9; % magnetic field (Tesla)

% initialize the parameters k, theta, kz, kx
par=zeros(6,1);
par(1)=0.5; % k, =sqrt(kx^2+kz^2), normalized by *c/omega_{p1}
par(2)=0; % theta, the angle between k and B0, normalized by *pi/180
par(3)=cos(par(2)*pi/180)*par(1); % kz, i.e., kpara*c/omega_{p1}
par(4)=sin(par(2)*pi/180)*par(1); % kx, i.e., kperp*c/omega_{p1}

% Choose which parameter(s) to scan, if ipa==ipb, do 1D scan; otherwise, 
% do 2D scan. % 1: k; 2: theta; 3: kz; 4: kx; 5: others.
% You can scan other parameters by modify line xxx in 'pdrk_si_kernel.m',
% e.g., to scan vs0(1), B0, betasz, etc
% Typical cases of (ipa,ipb): 
%   1. (1,1) scan k, fixed theta
%   2. (2,2) scan theta, fixed k
%   3. (1,2) scan 2D (k, theta)
%   4. (3,3) scan kz, fixed kx
%   5. (4,4) scan kx, fixed kz
%   6. (3,4) scan 2D (kz,kx)
%   7. (..,..) others, please modify 'pdrk_si_kernel.m'
ipa=1;
ipb=1;

iloga=0; % ilog=0, linear scale; =1, log scale, i.e., 10^(p1:dp:p2)
ilogb=0;
pa1=0.1; pa2=400; dpa=2; % 1st parameter a, depends on ipa
pb1=1; pb2=60; dpb=2.5; % 2nd parameter b, depends on ipb

% wether calculate polarizations (dEx,dEy,dEz,dBx,dBy,dBz) for select omega
iout=1; % =1, only (omega,gamma); =2, also calculate (E,B7


