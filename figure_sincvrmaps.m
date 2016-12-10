function figure_sincvrmaps(bids,subj)
	
	%read in sinusoid CVR maps (full data)
	sincvr_mag=read_avw([bids(subj).func(1).results(1).cvr_mag '_reg']);
	sincvr_magvar=read_avw([bids(subj).func(1).results(1).cvr_magvar '_reg']);
	sincvr_magrsd=read_avw([bids(subj).func(1).results(1).cvr_magrsd '_reg']);
	sincvr_pha=read_avw([bids(subj).func(1).results(1).cvr_pha 'n_reg']);
	sincvr_phavar=read_avw([bids(subj).func(1).results(1).cvr_phavar '_reg']);
	sincvr_zstat=read_avw([bids(subj).func(1).analysis(1).feat 'stats/zfstat1']);
	
	%read in sinusoid CVR maps (5 mins data)
	sin5cvr_mag=read_avw([bids(subj).func(1).results(3).cvr_mag '_reg']);
	sin5cvr_magvar=read_avw([bids(subj).func(1).results(3).cvr_magvar '_reg']);
	sin5cvr_magrsd=read_avw([bids(subj).func(1).results(3).cvr_magrsd '_reg']);
	sin5cvr_pha=read_avw([bids(subj).func(1).results(3).cvr_pha 'n_reg']);
	sin5cvr_phavar=read_avw([bids(subj).func(1).results(3).cvr_phavar '_reg']);
	sin5cvr_zstat=read_avw([bids(subj).func(1).analysis(3).feat 'stats/zfstat1']);
	
	%read in sinusoid CVR maps (3 mins data)
	sin3cvr_mag=read_avw([bids(subj).func(1).results(5).cvr_mag '_reg']);
	sin3cvr_magvar=read_avw([bids(subj).func(1).results(5).cvr_magvar '_reg']);
	sin3cvr_magrsd=read_avw([bids(subj).func(1).results(5).cvr_magrsd '_reg']);
	sin3cvr_pha=read_avw([bids(subj).func(1).results(5).cvr_pha 'n_reg']);
	sin3cvr_phavar=read_avw([bids(subj).func(1).results(5).cvr_phavar '_reg']);
	sin3cvr_zstat=read_avw([bids(subj).func(1).analysis(5).feat 'stats/zfstat1']);

	slices=[13];

	figure;
	smontage(sincvr_zstat(:,:,slices),1,1,[0 20]);
	colormap hot;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid zf-stat - 7 mins');
	
	figure;
	smontage(sin5cvr_zstat(:,:,slices),1,1,[0 20]);
	colormap hot;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid zf-stat - 5 mins');
	
	figure;
	smontage(sin3cvr_zstat(:,:,slices),1,1,[0 20]);
	colormap hot;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid zf-stat - 3 mins');

	figure;
	smontage(sincvr_mag(:,:,slices)./bids(subj).func(1).results(1).petco2delta.*100,1,1,[0 0.6]);
	colormap gray;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR magnitude - 7 mins');
	
	figure;
	smontage(sin5cvr_mag(:,:,slices)./bids(subj).func(1).results(3).petco2delta.*100,1,1,[0 0.6]);
	colormap gray;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR magnitude - 5 mins');
	
	figure;
	smontage(sin3cvr_mag(:,:,slices)./bids(subj).func(1).results(5).petco2delta.*100,1,1,[0 0.6]);
	colormap gray;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR magnitude - 3 mins');

	figure;
	smontage(sincvr_pha(:,:,slices),1,1,[-pi pi]);
	colormap hsv;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR phase - 7 mins');
	
	figure;
	smontage(sin5cvr_pha(:,:,slices),1,1,[-pi pi]);
	colormap hsv;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR phase - 5 mins');
	
	figure;
	smontage(sin3cvr_pha(:,:,slices),1,1,[-pi pi]);
	colormap hsv;
	set(gca,'yticklabel',[]);
	set(gca,'xticklabel',[]);
	set(gca,'ytick',[]);
	set(gca,'xtick',[]);
	title('Sinusoid CVR phase - 3 mins');
	
	%keyboard;