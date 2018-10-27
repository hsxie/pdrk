% 18-10-06 07:05 Hua-sheng XIE, huashengxie@gmail.com, FRI-ENN, China
% Ackn.: Richard Denton (Dartmouth), Xin Tao (USTC), Jin-song Zhao (PMO), 
% etc ...
% Use this file to output the results.

% iout=2;
% npl1=4;
npl1=npl;
if(iout==2)
  % calculate the polarizations for wws
  sp=1; nw=1;
  Pola=zeros(npa,npb,npl,20);
  wws2=wws;
  icalp=1; % tell kernel.m which can calculate polarization now
  for jpl=1:npl1
    % in eigs(), use wws(jpa,jpb,jpl) as intitial guess for each (pa,pb)
    run ../modules/pdrk_kernel; % to update
  end
end

% % To do: add group velocity, etc, 18-10-18 15:00
%%
close all;
if(iout==2)
h=figure('unit','normalized','Position',[0.01 0.25 0.5 0.6],...
  'DefaultAxesFontSize',15);
for jpl=1:npl1
    if(ipa==ipb) % plot 1D polarizations, to update. 2018-10-19 18:10
      subplot(331); hold on; box on;
      plot(pas,real(wws2(:,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]); %ylim([-2.5,1.0]);
      xlabel(strpa);ylabel('\omega_r/\omega_{c1}');
      
      subplot(332); hold on; box on;
      plot(pas,imag(wws2(:,1,jpl)),'-','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('\omega_i/\omega_{c1}');
      
      subplot(333); hold on; box on;
      plot(pas,Pola(:,1,jpl,7)./Pola(:,1,jpl,8),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('Energy E/Energy B');
      
      subplot(334); hold on; box on;
      plot(pas,real(Pola(:,1,jpl,2)./(1i*Pola(:,1,jpl,1))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      plot(pas,real(Pola(:,1,jpl,2)./(1i*Pola(:,1,jpl,1))),'-','Color',...
          pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_y/(iE_x)');
      
      subplot(335); hold on; box on;
      plot(pas,real(Pola(:,1,jpl,1)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola(:,1,jpl,1)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_x');
      legend('Re','Im'); legend('boxoff');
      
      subplot(336); hold on; box on;
      plot(pas,real(Pola(:,1,jpl,2)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola(:,1,jpl,2)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_y');
      
      subplot(337); hold on; box on;
      plot(pas,real(Pola(:,1,jpl,3)),'-','Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,imag(Pola(:,1,jpl,3)),'--','Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('E_z');
      
      subplot(338); hold on; box on;
      plot(pas,sqrt(c2)*real(Pola(:,1,jpl,4)),'-',...
          'Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,sqrt(c2)*imag(Pola(:,1,jpl,4)),'--',...
          'Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('c*B_x');
      
      subplot(339); hold on; box on;
      plot(pas,sqrt(c2)*real(Pola(:,1,jpl,6)),'-',...
          'Color',pltc(jpl,:),'linewidth',2);
      hold on;
      plot(pas,sqrt(c2)*imag(Pola(:,1,jpl,6)),'--',...
          'Color',pltc(jpl,:),'linewidth',2);
      xlim([min(pa),max(pa)]);
      xlabel(strpa);ylabel('c*B_z');
    else % plot 2D polarizations, to update
      
      wwjp=squeeze(wws2(:,:,jpl));
      Polajp=squeeze(Pola(:,:,jpl,:));
      
      subplot(221);
      surf(ppa,ppb,real(wwjp)); hold on; box on; set(gca,'BoxStyle','full');
      xlabel([strpa,',ilogx=',num2str(iloga)]);
      ylabel([strpb,',ilogy=',num2str(ilogb)]);axis tight;
      zlabel('\omega_r/\omega_{c1}');
      
      subplot(222);
      surf(ppa,ppb,imag(wwjp)); hold on; box on;set(gca,'BoxStyle','full');
      xlabel(strpa); ylabel(strpb); axis tight;
      zlabel('\omega_i/\omega_{c1}');
    
      EyoveriEx=Polajp(:,:,2)./(1i*Polajp(:,:,1));
      
      subplot(223);
      surf(ppa,ppb,real(EyoveriEx));
      hold on; box on;set(gca,'BoxStyle','full');
      xlabel(strpa); ylabel(strpb); axis tight;
      zlabel('Re[E_y/(iE_x)]');
      
      subplot(223);
      surf(ppa,ppb,imag(EyoveriEx));
      hold on; box on;set(gca,'BoxStyle','full');
      xlabel(strpa); ylabel(strpb); axis tight;
      zlabel('Im[E_y/(iE_x)]');
      
    end
end

print(gcf,'-dpng',[savepath,'fig_pdrk_',figstr,'_pola.png']);
savefig([savepath,'fig_pdrk_',figstr,'_pola.fig']);

end

filename=['out_pdrk_S=',num2str(S),'_J=',num2str(J),...
    '_N=',num2str(N),'_B0=',num2str(B0),'.mat'];

allvars=whos;
tosave=cellfun(@isempty, regexp({allvars.class}, ...
    '^matlab\.(ui|graphics)\.')); % save workspace exclude figure(s)
save([savepath,filename], allvars(tosave).name); 
