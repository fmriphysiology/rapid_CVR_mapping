function figure_compare_sin_sin_ev(bids)

	for subj=1:10
	
	%read in tor stat map
	zstat=read_avw([bids(subj).func(2).analysis(2).feat 'thresh_zstat1']);

	%read in sinusoid CVR maps (full data)
	[sincvr_mag dims]=read_avw([bids(subj).func(1).results(8).cvr_mag  ]);
	sincvr_mag=sincvr_mag./bids(subj).func(1).results(8).petco2delta.*100;
	sincvr_magvar=read_avw([bids(subj).func(1).results(8).cvr_magvar  ]);
	sincvr_fstat=read_avw([bids(subj).func(1).analysis(8).feat 'thresh_zfstat1']);

	%read in sinusoid CVR maps (6mins)
	sin6cvr_mag=read_avw([bids(subj).func(1).results(9).cvr_mag  ]);
	sin6cvr_mag=sin6cvr_mag./bids(subj).func(1).results(9).petco2delta.*100;
	sin6cvr_magvar=read_avw([bids(subj).func(1).results(9).cvr_magvar  ]);
	sin6cvr_fstat=read_avw([bids(subj).func(1).analysis(9).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (5mins)
	sin5cvr_mag=read_avw([bids(subj).func(1).results(10).cvr_mag  ]);
	sin5cvr_mag=sin5cvr_mag./bids(subj).func(1).results(10).petco2delta.*100;
	sin5cvr_magvar=read_avw([bids(subj).func(1).results(10).cvr_magvar  ]);
	sin5cvr_fstat=read_avw([bids(subj).func(1).analysis(10).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (4mins)
	sin4cvr_mag=read_avw([bids(subj).func(1).results(11).cvr_mag  ]);
	sin4cvr_mag=sin4cvr_mag./bids(subj).func(1).results(11).petco2delta.*100;
	sin4cvr_magvar=read_avw([bids(subj).func(1).results(11).cvr_magvar  ]);	
	sin4cvr_fstat=read_avw([bids(subj).func(1).analysis(11).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (3mins)
	sin3cvr_mag=read_avw([bids(subj).func(1).results(12).cvr_mag  ]);
	sin3cvr_mag=sin3cvr_mag./bids(subj).func(1).results(12).petco2delta.*100;
	sin3cvr_magvar=read_avw([bids(subj).func(1).results(12).cvr_magvar  ]);
	sin3cvr_fstat=read_avw([bids(subj).func(1).analysis(12).feat 'thresh_zfstat1']);

	%read in sinusoid CVR maps (2mins)
	sin2cvr_mag=read_avw([bids(subj).func(1).results(13).cvr_mag  ]);
	sin2cvr_mag=sin2cvr_mag./bids(subj).func(1).results(13).petco2delta.*100;
	sin2cvr_magvar=read_avw([bids(subj).func(1).results(13).cvr_magvar  ]);
	sin2cvr_fstat=read_avw([bids(subj).func(1).analysis(13).feat 'thresh_zfstat1']);
	
	%read in sinusoid CVR maps (1mins)
	sin1cvr_mag=read_avw([bids(subj).func(1).results(14).cvr_mag  ]);
	sin1cvr_mag=sin1cvr_mag./bids(subj).func(1).results(14).petco2delta.*100;
	sin1cvr_magvar=read_avw([bids(subj).func(1).results(14).cvr_magvar  ]);
	sin1cvr_fstat=read_avw([bids(subj).func(1).analysis(14).feat 'thresh_zfstat1']);
		
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
	sincvr_fstat(:,:,slicerm)=0;

	sin6cvr_mag(:,:,slicerm)=0;
	sin6cvr_magvar(:,:,slicerm)=0;
	sin6cvr_fstat(:,:,slicerm)=0;

	sin5cvr_mag(:,:,slicerm)=0;
	sin5cvr_magvar(:,:,slicerm)=0;
	sin5cvr_fstat(:,:,slicerm)=0;

	sin4cvr_mag(:,:,slicerm)=0;
	sin4cvr_magvar(:,:,slicerm)=0;
	sin4cvr_fstat(:,:,slicerm)=0;

	sin3cvr_mag(:,:,slicerm)=0;
	sin3cvr_magvar(:,:,slicerm)=0;
	sin3cvr_fstat(:,:,slicerm)=0;

	sin2cvr_mag(:,:,slicerm)=0;
	sin2cvr_magvar(:,:,slicerm)=0;
	sin2cvr_fstat(:,:,slicerm)=0;

	sin1cvr_mag(:,:,slicerm)=0;
	sin1cvr_magvar(:,:,slicerm)=0;
	sin1cvr_fstat(:,:,slicerm)=0;
			
	%rearrange matrices
	gms=gm(:);
	wms=wm(:);
	atlass=atlas(:);

	zstats=zstat(:);
	
	sincvr_mags=sincvr_mag(:);
	sincvr_magvars=sincvr_magvar(:);
	sincvr_fstats=sincvr_fstat(:);

	sin6cvr_mags=sin6cvr_mag(:);
	sin6cvr_magvars=sin6cvr_magvar(:);
	sin6cvr_fstats=sin6cvr_fstat(:);
	
	sin5cvr_mags=sin5cvr_mag(:);
	sin5cvr_magvars=sin5cvr_magvar(:);
	sin5cvr_fstats=sin5cvr_fstat(:);
	
	sin4cvr_mags=sin4cvr_mag(:);
	sin4cvr_magvars=sin4cvr_magvar(:);
	sin4cvr_fstats=sin4cvr_fstat(:);
		
	sin3cvr_mags=sin3cvr_mag(:);
	sin3cvr_magvars=sin3cvr_magvar(:);
	sin3cvr_fstats=sin3cvr_fstat(:);
	
	sin2cvr_mags=sin2cvr_mag(:);
	sin2cvr_magvars=sin2cvr_magvar(:);
	sin2cvr_fstats=sin2cvr_fstat(:);
	
	sin1cvr_mags=sin1cvr_mag(:);
	sin1cvr_magvars=sin1cvr_magvar(:);
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
	
%	figure; %sd CVR mag & pha normalised to 7mins of data
%	errorbar((7:-1:1),mean(allsd(1:9,:)./mean(allsd(1:9,1))),std(allsd(1:9,:)./mean(allsd(1:9,1))),'o-');
%	hold on;
%	errorbar((7:-1:1),mean(allphasd(1:9,:)./mean(allphasd(1:9,1))),std(allphasd(1:9,:)./mean(allphasd(1:9,1))),'o-');
%	axis square;
%	ylim([0.6 2.4]);
%	set(gca,'ytick',[0.6:0.4:2.2]);
%	title('Standard deviation of CVR Mag/Phase normalised to 7 minutes');
	
	