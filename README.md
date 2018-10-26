pdrk
====
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

Next version will be more user-friendly and support loss-cone drift bi-Maxwellian distribution.

2018-10-03 13:14

----------------------------------------------------------------

Plasma Dispersion Relation Solver , Kinetic Version

Supplementary material for the MATLAB numerical routines used in PDRK paper.

Author: H. S. Xie & Y. Xiao, huashengxie@gmail.com, IFTS-ZJU

Title: PDRK: A General Kinetic Dispersion Relation Solver for Magnetized Plasma

Journal: XXX (or see arxiv: http://arxiv.org/abs/1410.2678)

Files:

./Cal_J_pole_bjcj.m    -- Calculate J-pole (Pade approximation) coefficient
                        bj and cj for Z_J(zeta)=sum(bj/(zeta-cj))
						
./pdrk_es1d.m          -- PDRK solver, ES1D
./pdrk_es1d.in         -- input file

./pdrk_es3d.m          -- PDRK solver, ES3D
./pdrk_es3d.in         -- input file

./pdrk_em3d.m          -- PDRK solver, EM3D
./pdrk_em3d.in         -- input file


2014-10-10 13:21
