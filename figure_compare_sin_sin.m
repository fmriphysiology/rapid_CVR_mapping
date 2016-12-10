function figure_compare_sin_sin(bids)

	for subj=1:10
	
	%read in tor stat map
	zstat=read_avw([bids(subj).func(2).analysis(2).feat 'thresh_zstat1']);

	%read in sinusoid CVR maps (full data)
	[sincvr_mag dims]=read_avw([bids(subj).func(1).results(1).cvr_mag '_reg']);
	sincvr_mag=sincvr_mag./bids(subj).func(1).results(1).petco2delta.*100;
	sincvr_magvar=read_avw([bids(subj).func(1).results(1).cvr_magvar '_reg']);
	sincvr_pha=read_avw([bids(subj).func(1).results(1).cvr_pha 'n_reg']);
	sincvr_phavar=read_avw([bids(subj).func(1).results(1).cvr_phavar '_reg']);
	sincvr_fstat=read_avw([bids(subj).func(1).analysis(1).feat 'thresh_zfstat1']);

	%read in sinusoid CVR maps (6mins)
	sin6cvr_mag=read_avw([bids(subj).func(1).results(2).cvr_mag '_reg']);
	sin6cvr_mag=sin6cvr_mag./bids(subj).func(1).results(2).petco2delta.*100;
	sin6cvr_magvar=read_avw([bids(subj).func(1).results(2).cvr_magvar '_reg']);
	sin6cvr_pha=read_avw([bids(subj).func(1).results(2).cvr_pha 'n_reg']);
	sin6cvr_phavar=read_avw([bids(subj).func(1).results(2).cvr_phavar '_reg']);
	sin6cvr_fstat=read_avw([bids(subj).func(1).analysis(2).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (5mins)
	sin5cvr_mag=read_avw([bids(subj).func(1).results(3).cvr_mag '_reg']);
	sin5cvr_mag=sin5cvr_mag./bids(subj).func(1).results(3).petco2delta.*100;
	sin5cvr_magvar=read_avw([bids(subj).func(1).results(3).cvr_magvar '_reg']);
	sin5cvr_pha=read_avw([bids(subj).func(1).results(3).cvr_pha 'n_reg']);
	sin5cvr_phavar=read_avw([bids(subj).func(1).results(3).cvr_phavar '_reg']);
	sin5cvr_fstat=read_avw([bids(subj).func(1).analysis(3).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (4mins)
	sin4cvr_mag=read_avw([bids(subj).func(1).results(4).cvr_mag '_reg']);
	sin4cvr_mag=sin4cvr_mag./bids(subj).func(1).results(4).petco2delta.*100;
	sin4cvr_magvar=read_avw([bids(subj).func(1).results(4).cvr_magvar '_reg']);
	sin4cvr_pha=read_avw([bids(subj).func(1).results(4).cvr_pha 'n_reg']);
	sin4cvr_phavar=read_avw([bids(subj).func(1).results(4).cvr_phavar '_reg']);	
	sin4cvr_fstat=read_avw([bids(subj).func(1).analysis(4).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (3mins)
	sin3cvr_mag=read_avw([bids(subj).func(1).results(5).cvr_mag '_reg']);
	sin3cvr_mag=sin3cvr_mag./bids(subj).func(1).results(5).petco2delta.*100;
	sin3cvr_magvar=read_avw([bids(subj).func(1).results(5).cvr_magvar '_reg']);
	sin3cvr_pha=read_avw([bids(subj).func(1).results(5).cvr_pha 'n_reg']);
	sin3cvr_phavar=read_avw([bids(subj).func(1).results(5).cvr_phavar '_reg']);
	sin3cvr_fstat=read_avw([bids(subj).func(1).analysis(5).feat 'thresh_zfstat1']);

	%read in sinusoid CVR maps (2mins)
	sin2cvr_mag=read_avw([bids(subj).func(1).results(6).cvr_mag '_reg']);
	sin2cvr_mag=sin2cvr_mag./bids(subj).func(1).results(6).petco2delta.*100;
	sin2cvr_magvar=read_avw([bids(subj).func(1).results(6).cvr_magvar '_reg']);
	sin2cvr_pha=read_avw([bids(subj).func(1).results(6).cvr_pha 'n_reg']);
	sin2cvr_phavar=read_avw([bids(subj).func(1).results(6).cvr_phavar '_reg']);
	sin2cvr_fstat=read_avw([bids(subj).func(1).analysis(6).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (1mins)
	sin1cvr_mag=read_avw([bids(subj).func(1).results(7).cvr_mag '_reg']);
	sin1cvr_mag=sin1cvr_mag./bids(subj).func(1).results(7).petco2delta.*100;
	sin1cvr_magvar=read_avw([bids(subj).func(1).results(7).cvr_magvar '_reg']);
	sin1cvr_pha=read_avw([bids(subj).func(1).results(7).cvr_pha 'n_reg']);
	sin1cvr_phavar=read_avw([bids(subj).func(1).results(7).cvr_phavar '_reg']);
	sin1cvr_fstat=read_avw([bids(subj).func(1).analysis(7).feat 'thresh_zfstat1']);
		
	%read in tissue masks
	gm=read_avw([bids(subj).anat.gm '_reg']);
	wm=read_avw([bids(subj).anat.wm '_reg']);
	
	s=regexp(bids(subj).func(1).results(1).cvr_mag,'/');
	out=[bids(subj).func(1).results(1).cvr_mag(1:s(end)) 'atlas_reg'];
	atlas=read_avw(out);
	
	%discard top and bottom slices (lost to realignment)
	slicerm=[1 dims(3)];
		
	gm(:,:,slicerm)=0;
	wm(:,:,slicerm)=0;
	atlas(:,:,slicerm)=0;

	zstat(:,:,slicerm)=0;
	
	sincvr_mag(:,:,slicerm)=0;
	sincvr_magvar(:,:,slicerm)=0;
	sincvr_pha(:,:,slicerm)=0;
	sincvr_phavar(:,:,slicerm)=0;
	sincvr_fstat(:,:,slicerm)=0;

	sin6cvr_mag(:,:,slicerm)=0;
	sin6cvr_magvar(:,:,slicerm)=0;
	sin6cvr_pha(:,:,slicerm)=0;
	sin6cvr_phavar(:,:,slicerm)=0;
	sin6cvr_fstat(:,:,slicerm)=0;

	sin5cvr_mag(:,:,slicerm)=0;
	sin5cvr_magvar(:,:,slicerm)=0;
	sin5cvr_pha(:,:,slicerm)=0;
	sin5cvr_phavar(:,:,slicerm)=0;
	sin5cvr_fstat(:,:,slicerm)=0;

	sin4cvr_mag(:,:,slicerm)=0;
	sin4cvr_magvar(:,:,slicerm)=0;
	sin4cvr_pha(:,:,slicerm)=0;
	sin4cvr_phavar(:,:,slicerm)=0;
	sin4cvr_fstat(:,:,slicerm)=0;

	sin3cvr_mag(:,:,slicerm)=0;
	sin3cvr_magvar(:,:,slicerm)=0;
	sin3cvr_pha(:,:,slicerm)=0;
	sin3cvr_phavar(:,:,slicerm)=0;
	sin3cvr_fstat(:,:,slicerm)=0;

	sin2cvr_mag(:,:,slicerm)=0;
	sin2cvr_magvar(:,:,slicerm)=0;
	sin2cvr_pha(:,:,slicerm)=0;
	sin2cvr_phavar(:,:,slicerm)=0;
	sin2cvr_fstat(:,:,slicerm)=0;

	sin1cvr_mag(:,:,slicerm)=0;
	sin1cvr_magvar(:,:,slicerm)=0;
	sin1cvr_pha(:,:,slicerm)=0;
	sin1cvr_phavar(:,:,slicerm)=0;
	sin1cvr_fstat(:,:,slicerm)=0;
			
	%rearrange matrices
	gms=gm(:);
	wms=wm(:);
	atlass=atlas(:);

	zstats=zstat(:);
	
	sincvr_mags=sincvr_mag(:);
	sincvr_magvars=sincvr_magvar(:);
	sincvr_phas=sincvr_pha(:);
	sincvr_phavars=sincvr_phavar(:);
	sincvr_fstats=sincvr_fstat(:);

	sin6cvr_mags=sin6cvr_mag(:);
	sin6cvr_magvars=sin6cvr_magvar(:);
	sin6cvr_phas=sin6cvr_pha(:);
	sin6cvr_phavars=sin6cvr_phavar(:);
	sin6cvr_fstats=sin6cvr_fstat(:);
	
	sin5cvr_mags=sin5cvr_mag(:);
	sin5cvr_magvars=sin5cvr_magvar(:);
	sin5cvr_phas=sin5cvr_pha(:);
	sin5cvr_phavars=sin5cvr_phavar(:);
	sin5cvr_fstats=sin5cvr_fstat(:);
	
	sin4cvr_mags=sin4cvr_mag(:);
	sin4cvr_magvars=sin4cvr_magvar(:);
	sin4cvr_phas=sin4cvr_pha(:);
	sin4cvr_phavars=sin4cvr_phavar(:);	
	sin4cvr_fstats=sin4cvr_fstat(:);
		
	sin3cvr_mags=sin3cvr_mag(:);
	sin3cvr_magvars=sin3cvr_magvar(:);
	sin3cvr_phas=sin3cvr_pha(:);
	sin3cvr_phavars=sin3cvr_phavar(:);	
	sin3cvr_fstats=sin3cvr_fstat(:);
	
	sin2cvr_mags=sin2cvr_mag(:);
	sin2cvr_magvars=sin2cvr_magvar(:);
	sin2cvr_phas=sin2cvr_pha(:);
	sin2cvr_phavars=sin2cvr_phavar(:);	
	sin2cvr_fstats=sin2cvr_fstat(:);
	
	sin1cvr_mags=sin1cvr_mag(:);
	sin1cvr_magvars=sin1cvr_magvar(:);
	sin1cvr_phas=sin1cvr_pha(:);
	sin1cvr_phavars=sin1cvr_phavar(:);	
	sin1cvr_fstats=sin1cvr_fstat(:);
	
	%define mask
	mask=(zstats>0).*(gms>0.5); %significant response to toronto stimulus in grey matter
	
	sincvr_magavg(subj,:)=mean(sincvr_mags(mask>0));
	sin6cvr_magavg(subj,:)=mean(sin6cvr_mags(mask>0));
	sin5cvr_magavg(subj,:)=mean(sin5cvr_mags(mask>0));
	sin4cvr_magavg(subj,:)=mean(sin4cvr_mags(mask>0));
	sin3cvr_magavg(subj,:)=mean(sin3cvr_mags(mask>0));
	sin2cvr_magavg(subj,:)=mean(sin2cvr_mags(mask>0));
	sin1cvr_magavg(subj,:)=mean(sin1cvr_mags(mask>0));
	
	sincvr_magsd(subj,:)=std(sincvr_mags(mask>0));
	sin6cvr_magsd(subj,:)=std(sin6cvr_mags(mask>0));
	sin5cvr_magsd(subj,:)=std(sin5cvr_mags(mask>0));
	sin4cvr_magsd(subj,:)=std(sin4cvr_mags(mask>0));
	sin3cvr_magsd(subj,:)=std(sin3cvr_mags(mask>0));
	sin2cvr_magsd(subj,:)=std(sin2cvr_mags(mask>0));
	sin1cvr_magsd(subj,:)=std(sin1cvr_mags(mask>0));

	sincvr_phaavg(subj,:)=mean(sincvr_phas(mask>0));
	sin6cvr_phaavg(subj,:)=mean(sin6cvr_phas(mask>0));
	sin5cvr_phaavg(subj,:)=mean(sin5cvr_phas(mask>0));
	sin4cvr_phaavg(subj,:)=mean(sin4cvr_phas(mask>0));
	sin3cvr_phaavg(subj,:)=mean(sin3cvr_phas(mask>0));
	sin2cvr_phaavg(subj,:)=mean(sin2cvr_phas(mask>0));
	sin1cvr_phaavg(subj,:)=mean(sin1cvr_phas(mask>0));
	
	sincvr_phasd(subj,:)=std(sincvr_phas(mask>0));
	sin6cvr_phasd(subj,:)=std(sin6cvr_phas(mask>0));
	sin5cvr_phasd(subj,:)=std(sin5cvr_phas(mask>0));
	sin4cvr_phasd(subj,:)=std(sin4cvr_phas(mask>0));
	sin3cvr_phasd(subj,:)=std(sin3cvr_phas(mask>0));
	sin2cvr_phasd(subj,:)=std(sin2cvr_phas(mask>0));
	sin1cvr_phasd(subj,:)=std(sin1cvr_phas(mask>0));

	sincvr_statnum(subj,:)=sum(sincvr_fstats>0);
	sin6cvr_statnum(subj,:)=sum(sin6cvr_fstats>0);
	sin5cvr_statnum(subj,:)=sum(sin5cvr_fstats>0);
	sin4cvr_statnum(subj,:)=sum(sin4cvr_fstats>0);
	sin3cvr_statnum(subj,:)=sum(sin3cvr_fstats>0);
	sin2cvr_statnum(subj,:)=sum(sin2cvr_fstats>0);
	sin1cvr_statnum(subj,:)=sum(sin1cvr_fstats>0);
	
	end
	
	allavg=[sincvr_magavg sin6cvr_magavg sin5cvr_magavg sin4cvr_magavg sin3cvr_magavg sin2cvr_magavg sin1cvr_magavg];
	allsd=[sincvr_magsd sin6cvr_magsd sin5cvr_magsd sin4cvr_magsd sin3cvr_magsd sin2cvr_magsd sin1cvr_magsd];

	allphaavg=[sincvr_phaavg sin6cvr_phaavg sin5cvr_phaavg sin4cvr_phaavg sin3cvr_phaavg sin2cvr_phaavg sin1cvr_phaavg];
	allphasd=[sincvr_phasd sin6cvr_phasd sin5cvr_phasd sin4cvr_phasd sin3cvr_phasd sin2cvr_phasd sin1cvr_phasd];
	
	allstatnum=[sincvr_statnum sin6cvr_statnum sin5cvr_statnum sin4cvr_statnum sin3cvr_statnum sin2cvr_statnum sin1cvr_statnum];
	
	%figures - ignore subj 10
	
	figure; %mean CVR mag vs time
	errorbar((7:-1:1),mean(allavg(1:9,:)),std(allavg(1:9,:)),'o-')
	axis square;
	ylim([0.3 0.8]);
	set(gca,'ytick',[0.3:0.1:0.8])
	title('Mean CVR Magnitude as a function of scan duration');
	
	[p table stats]=anova2(allavg(1:9,:));
	if p(1)<0.05
		figure;
		comparison=multcompare(stats);
		s=comparison(find(comparison(find(comparison(:,1)==1),end)<0.05),2);
		fprintf('CVR Magnitude: Full data is significantly different to the cycles:');
		for k=1:length(s)
			fprintf([num2str(8-s(k)) ' '])
		end
		fprintf('\n');
	else
		disp('CVR Magnitude: Groups are not significantly different');
	end
	
	figure; %sd CVR mag vs time
	errorbar((7:-1:1),mean(allsd(1:9,:)),std(allsd(1:9,:)),'o-')
	axis square;
	ylim([0.3 0.9]);
	set(gca,'ytick',[0.3:0.1:0.9])	 
	title('Standard deviation of CVR Magnitude as a function of scan duration');

	[p table stats]=anova2(allsd(1:9,:));
	if p(1)<0.05
		figure;
		comparison=multcompare(stats);
		s=comparison(find(comparison(find(comparison(:,1)==1),end)<0.05),2);
		fprintf('CVR Magnitude SD: Full data is significantly different to the cycles:');
		for k=1:length(s)
			fprintf([num2str(8-s(k)) ' '])
		end
		fprintf('\n');
	else
		disp('CVR Magnitude SD: Groups are not significantly different');
	end

	figure; %mean CVR pha vs time
	errorbar((7:-1:1),mean(allphaavg(1:9,:)),std(allphaavg(1:9,:)),'o-')
	axis square;
	ylim([-0.3 0.3]);
	set(gca,'ytick',[-0.2:0.1:0.3])
	title('Mean CVR Phase as a function of scan duration');

	[p table stats]=anova2(allphaavg(1:9,:));
	if p(1)<0.05
		figure;
		comparison=multcompare(stats);
		s=comparison(find(comparison(find(comparison(:,1)==1),end)<0.05),2);
		fprintf('CVR Phase: Full data is significantly different to the cycles:');
		for k=1:length(s)
			fprintf([num2str(8-s(k)) ' '])
		end
		fprintf('\n');
	else
		disp('CVR Phase: Groups are not significantly different');
	end
	
	figure; %sd CVR pha vs time
	errorbar((7:-1:1),mean(allphasd(1:9,:)),std(allphasd(1:9,:)),'o-')
	axis square;
	ylim([0.1 0.6]);
	set(gca,'ytick',[0.1:0.1:0.6])
	title('Standard deviation of CVR Phase as a function of scan duration');
	
	[p table stats]=anova2(allphasd(1:9,:));
	if p(1)<0.05
		figure;
		comparison=multcompare(stats);
		s=comparison(find(comparison(find(comparison(:,1)==1),end)<0.05),2);
		fprintf('CVR Magnitude SD: Full data is significantly different to the cycles:');
		for k=1:length(s)
			fprintf([num2str(8-s(k)) ' '])
		end
		fprintf('\n');
	else
		disp('CVR Magnitude SD: Groups are not significantly different');
	end
	
	figure; %mean number of voxels vs time
	errorbar((7:-1:1),mean(allstatnum(1:9,:)),std(allstatnum(1:9,:)),'o-')
	axis square;
	ylim([0 2e4]);
	set(gca,'ytick',[0:0.5e4:2e4]);
	title('Number of voxels above threshold as a function of scan duration');

	[p table stats]=anova2(allstatnum(1:9,:));
	if p(1)<0.05
		figure;
		comparison=multcompare(stats);
		s=comparison(find(comparison(find(comparison(:,1)==1),end)<0.05),2);
		fprintf('CVR Statistical Num Voxels: Full data is significantly different to the cycles:');
		for k=1:length(s)
			fprintf([num2str(8-s(k)) ' '])
		end
		fprintf('\n');
	else
		disp('CVR Magnitude SD: Groups are not significantly different');
	end
	
	figure; %sd CVR mag & pha normalised to 7mins of data
	errorbar((7:-1:1),mean(allsd(1:9,:)./mean(allsd(1:9,1))),std(allsd(1:9,:)./mean(allsd(1:9,1))),'o-');
	hold on;
	errorbar((7:-1:1),mean(allphasd(1:9,:)./mean(allphasd(1:9,1))),std(allphasd(1:9,:)./mean(allphasd(1:9,1))),'o-');
	axis square;
	ylim([0.6 2.4]);
	set(gca,'ytick',[0.6:0.4:2.2]);
	title('Standard deviation of CVR Mag/Phase normalised to 7 minutes');
	
	