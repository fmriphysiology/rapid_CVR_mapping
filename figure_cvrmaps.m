function figure_cvrmaps(bids,subj)

	%read in toronto CVR maps 
	[torcvr_mag dims scales]=read_avw(bids(subj).func(2).results(1).cvr_mag);
	torcvr_magvar=read_avw(bids(subj).func(2).results(1).cvr_magvar);
	torcvr_magrsd=read_avw(bids(subj).func(2).results(1).cvr_magrsd);
	torcvr_pha=read_avw(bids(subj).func(2).results(1).cvr_pha);
	torcvr_zstat=read_avw([bids(subj).func(2).analysis(2).feat 'stats/zstat1']);
	
	%read in sinusoid CVR maps (full data)
	sincvr_mag=read_avw([bids(subj).func(1).results(1).cvr_mag '_reg']);
	sincvr_magvar=read_avw([bids(subj).func(1).results(1).cvr_magvar '_reg']);
	sincvr_magrsd=read_avw([bids(subj).func(1).results(1).cvr_magrsd '_reg']);
	sincvr_pha=read_avw([bids(subj).func(1).results(1).cvr_pha 'n_reg']);
	sincvr_phavar=read_avw([bids(subj).func(1).results(1).cvr_phavar '_reg']);
	sincvr_zstat=read_avw([bids(subj).func(1).analysis(1).feat 'stats/zfstat1']);

	slices=[13:17];

	figure;
	smontage(torcvr_zstat(:,:,slices),1,5,[0 20]);
	colormap hot;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Toronto z-stat');

	figure;
	smontage(sincvr_zstat(:,:,slices),1,5,[0 20]);
	colormap hot;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid zf-stat');

	figure;	
	smontage(abs(torcvr_mag(:,:,slices))./bids(subj).func(2).results(1).petco2delta.*100,1,5,[0 0.6]);
	colormap gray;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	h=colorbar('fontsize',16);
	set(h,'ticks',[0 3 6]);
	title('Toronto CVR Magnitude');

	figure;
	smontage(sincvr_mag(:,:,slices)./bids(subj).func(1).results(1).petco2delta.*100,1,5,[0 0.6]);
	colormap gray;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	h=colorbar('fontsize',16);
	set(h,'ticks',[0 3 6]);
	title('Sinusoid CVR Magnitude');

	figure;	
	smontage(torcvr_magrsd(:,:,slices).*100,1,5,[0 40]);
	colormap jet;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Toronto CVR Magnitude Relative Standard Deviation');

	figure;	
	smontage(sincvr_magrsd(:,:,slices).*100,1,5,[0 40]);
	colormap jet;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR Magnitude Relative Standard Deviation');

	figure;	
	smontage(torcvr_pha(:,:,slices),1,5,[-pi pi]);
	colormap hsv;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	h=colorbar('fontsize',16);
	set(h,'ticks',[-pi 0 pi]);
	set(h,'ticklabels',{'-\pi' '0' '\pi'});
	title('Toronto CVR Phase');
	
	figure;	
	smontage(sincvr_pha(:,:,slices),1,5,[-pi pi]);
	colormap hsv;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	h=colorbar('fontsize',16);
	set(h,'ticks',[-pi 0 pi]);
	set(h,'ticklabels',{'-\pi' '0' '\pi'});
	title('Sinusoid CVR Phase');
	
	figure;	
	smontage(sqrt(sincvr_phavar(:,:,slices)),1,5,[0 0.4]);
	colormap jet;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR Phase Standard Deviation');
	