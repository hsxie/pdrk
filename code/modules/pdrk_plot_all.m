% 18-10-05 08:00 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO),
% etc ...
% This file plot all the solution of pdrk-em3d results, which can be used
% to select the initial data for separating different dispersion surface
close all;

h=figure('unit','normalized','Position',[0.01 0.45 0.5 0.3],...
  'DefaultAxesFontSize',15);

% normalized omega and k
wwn=ww/wcs1;
kkn=kk*cwp;
kxxn=kxx*cwp;
kzzn=kzz*cwp;

if(ipa==ipb) % 1D plot

  %kkk=reshape(repmat(kkn,1,nw0),1,[]);
  papp=reshape(repmat(ppa,1,nw0),1,[]);
  www=reshape(wwn,1,[]);
  npp=length(papp);

  subplot(121);
  if(iloga==0)
    plot(papp,real(www),'g+','LineWidth',2); hold on;
  else
    semilogx(10.^papp,real(www),'g+','LineWidth',2); hold on;
  end  
  xlabel([strpa,', runtime=',num2str(runtime),'s']); 
  ylabel(['\omega_r/\omega_{c1}, npa=',num2str(npa),',npb=',num2str(npb)]);
  title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]);
  xlim([min(pa),max(pa)]); box on; %ylim([-2.5,2.5]);
  
  subplot(122);
  if(iloga==0)
    plot(papp,imag(www),'g+','LineWidth',2); hold on;
  else
    semilogx(10.^papp,imag(www),'g+','LineWidth',2); hold on;
  end
  xlabel([strpa,', (S=',num2str(S),',N=',num2str(N),',J=',num2str(J),')']); 
  ylabel('\omega_i/\omega_{c1}');
  title(['(b) v_A/c=',num2str(vA/sqrt(c2),2),', ',strpb,'=',...
    num2str(par(ipbtmp))]);
  xlim([min(pa),max(pa)]); box on; %ylim([-1.0,0.1]);
  
else % 2D plot
  for jp=1:25 % change here to plot more surfaces
    subplot(121);
    wwjp=squeeze(wwn(:,:,jp));
    surf(ppa,ppb,real(wwjp)); hold on; box on;
    xlabel([strpa,',ilogx=',num2str(iloga)]);
    ylabel([strpb,',ilogy=',num2str(ilogb)]);
    zlabel(['\omega_r/\omega_{c1},npa=',num2str(npa),',npb=',num2str(npb)]);
    title(['(a) \beta_{||}=',num2str(betasz,3),...
      ', \beta_\perp=',num2str(betasp,3)]);
    subplot(122);
    surf(ppa,ppb,imag(wwjp)); hold on; box on;
    xlabel(strpa); ylabel(strpb);
    zlabel(['\omega_i/\omega_{c1},N=',num2str(N),',J=',num2str(J)]);
    title(['(b) runtime=',num2str(runtime),'s']);
    %%
    % zoom in the figure to find dispersion surface data for plot_select.m
    %zlim([-0.5e0,0.1]);
    %%
  end
  
end

figstr=['S=',num2str(S),'_J=',num2str(J),'_N=',num2str(N),...
    '_npa=',num2str(npa),'_npb=',num2str(npb)];
print(gcf,'-dpng',[savepath,'fig_pdrk_',figstr,'_all.png']);
savefig([savepath,'fig_pdrk_',figstr,'_all.fig']);
