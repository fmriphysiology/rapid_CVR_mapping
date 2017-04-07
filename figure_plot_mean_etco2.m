function figure_plot_mean_etco2(bids)

	for subj=1:9
		
		[gm(:,:,subj) wm(:,:,subj) t1 petco2(:,:,subj) t2]=plot_bold_etco2(bids,subj);
		
	end
	
	figure;
	subplot(311);
	plot(t2./60,mean(squeeze(petco2(:,1,:)),2),'k');
	hold on;
	plot(t2./60,mean(squeeze(petco2(:,1,:)),2)+std(squeeze(petco2(:,1,:)),[],2),'k--');
	plot(t2./60,mean(squeeze(petco2(:,1,:)),2)-std(squeeze(petco2(:,1,:)),[],2),'k--');
	xlim([-1 8]);
	ylim([30 55]);
	grid on;
	title('Sinusoid: Mean end-tidal CO_2');
	
	subplot(312);
	plot(t1./60,mean(squeeze(gm(:,1,:)),2),'k');
	hold on;
	plot(t1./60,mean(squeeze(gm(:,1,:)),2)+std(squeeze(gm(:,1,:)),[],2),'k--');
	plot(t1./60,mean(squeeze(gm(:,1,:)),2)-std(squeeze(gm(:,1,:)),[],2),'k--');
	xlim([-1 8]);
	ylim([-4 4]);
	grid on;
	title('Sinusoid: Mean grey matter BOLD signal');
	
	subplot(313);
	plot(t1./60,mean(squeeze(wm(:,1,:)),2),'k');
	hold on;
	plot(t1./60,mean(squeeze(wm(:,1,:)),2)+std(squeeze(wm(:,1,:)),[],2),'k--');
	plot(t1./60,mean(squeeze(wm(:,1,:)),2)-std(squeeze(wm(:,1,:)),[],2),'k--');
	xlim([-1 8]);
	ylim([-2 2]);	
	grid on;
	title('Sinusoid: Mean white matter BOLD signal');
	
	figure;
	subplot(311);
	plot(t2./60,mean(squeeze(petco2(:,2,:)),2),'k');
	hold on;
	plot(t2./60,mean(squeeze(petco2(:,2,:)),2)+std(squeeze(petco2(:,2,:)),[],2),'k--');
	plot(t2./60,mean(squeeze(petco2(:,2,:)),2)-std(squeeze(petco2(:,2,:)),[],2),'k--');
	xlim([-1 8]);
	ylim([30 55]);
	grid on;
	title('Toronto: Mean end-tidal CO_2');
	
	subplot(312);
	plot(t1./60,mean(squeeze(gm(:,2,:)),2),'k');
	hold on;
	plot(t1./60,mean(squeeze(gm(:,2,:)),2)+std(squeeze(gm(:,2,:)),[],2),'k--');
	plot(t1./60,mean(squeeze(gm(:,2,:)),2)-std(squeeze(gm(:,2,:)),[],2),'k--');
	xlim([-1 8]);
	ylim([-4 4]);
	grid on;
	title('Toronto: Mean grey matter BOLD signal');
		
	subplot(313);
	plot(t1./60,mean(squeeze(wm(:,2,:)),2),'k');
	hold on;
	plot(t1./60,mean(squeeze(wm(:,2,:)),2)+std(squeeze(wm(:,2,:)),[],2),'k--');
	plot(t1./60,mean(squeeze(wm(:,2,:)),2)-std(squeeze(wm(:,2,:)),[],2),'k--');
	xlim([-1 8]);
	ylim([-2 2]);	
	grid on;
	title('Toronto: Mean white matter BOLD signal');
		
	%keyboard;
		