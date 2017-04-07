function figure_compare_tor_sin_ev(bids)

	%subjects to include in analysis
	subj_inc=1:9; 
	%example subject for displaying plot
	exsubj=1;

	options=optimset('fminsearch');

	for subj=1:10

	%read in toronto CVR maps 
	[torcvr_mag dims scales]=read_avw(bids(subj).func(2).results(1).cvr_mag);
	torcvr_magvar=read_avw(bids(subj).func(2).results(1).cvr_magvar);
	torcvr_pha=read_avw(bids(subj).func(2).results(1).cvr_pha);

	%read in sinusoid CVR maps (full data)
	sincvr_mag=read_avw([bids(subj).func(1).results(8).cvr_mag]);
	sincvr_magvar=read_avw([bids(subj).func(1).results(8).cvr_magvar]);
	sincvr_pha=read_avw([bids(subj).func(1).results(1).cvr_pha '_reg']);
	sincvr_phavar=read_avw([bids(subj).func(1).results(1).cvr_phavar '_reg']);
	
	%read in tissue masks
	gm=read_avw([bids(subj).anat.gm '_reg']);
	
	%discard top and bottom slices (lost to realignment)
	slicerm=[1 dims(3)];
		
	gm(:,:,slicerm)=0;
	
	torcvr_mag(:,:,slicerm)=0;
	torcvr_magvar(:,:,slicerm)=0;
	torcvr_pha(:,:,slicerm)=0;
	
	sincvr_mag(:,:,slicerm)=0;
	sincvr_magvar(:,:,slicerm)=0;
	sincvr_pha(:,:,slicerm)=0;
	sincvr_phavar(:,:,slicerm)=0;
	
	%rearrange matrices
	gms=gm(:);
	
	torcvr_mags=torcvr_mag(:);
	torcvr_magvars=torcvr_magvar(:);
	torcvr_phas=torcvr_pha(:);
	
	sincvr_mags=sincvr_mag(:);
	sincvr_magvars=sincvr_magvar(:);
	sincvr_phas=sincvr_pha(:);
	sincvr_phavars=sincvr_phavar(:);
	
	torcvr_cvs=sqrt(torcvr_magvars)./torcvr_mags;
	sincvr_cvs=sqrt(sincvr_magvars)./sincvr_mags;
	
	%define mask
	mask=(sincvr_mags~=0).*(torcvr_mags~=0).*(gms>0.5);
	
	%normalise magnitude by petco2 change
	torcvr_mags=abs(torcvr_mags)./bids(subj).func(2).results(1).petco2delta;
	sincvr_mags=sincvr_mags./bids(subj).func(1).results(1).petco2delta;

	%normalise phase by mean ROI cvr phase value
	sincvr_phasn=sincvr_phas;
	sincvr_phasn(sincvr_phasn~=0)=sincvr_phasn(sincvr_phasn~=0)-mean(sincvr_phasn(mask>0));
	sincvr_phasn(sincvr_phasn>pi)=sincvr_phasn(sincvr_phasn>pi)-2*pi;
	sincvr_phasn(sincvr_phasn<-pi)=sincvr_phasn(sincvr_phasn<-pi)+2*pi;
	
	torcvr_phasn=torcvr_phas;
	torcvr_phasn(torcvr_phasn~=0)=torcvr_phasn(torcvr_phasn~=0)-mean(torcvr_phasn(mask>0));
	torcvr_phasn(torcvr_phasn>pi)=torcvr_phasn(torcvr_phasn>pi)-2*pi;
	torcvr_phasn(torcvr_phasn<-pi)=torcvr_phasn(torcvr_phasn<-pi)+2*pi;
	
	Xpha=[sin(torcvr_phasn(mask>0)) ones(size(sin(torcvr_phasn(mask>0))))];
	a_pha(:,subj)=Xpha\sincvr_phasn(mask>0);
	
	Xpha=[sin(sincvr_phasn(mask>0)) ones(size(sin(torcvr_phasn(mask>0))))];
	a_pha2(:,subj)=Xpha\torcvr_phasn(mask>0);
	
	[corr_pha(:,:,subj) sig_pha(:,:,subj)]=corrcoef(sin(torcvr_phasn(mask>0)),sin(sincvr_phasn(mask>0)));
	
	diffmap=torcvr_phasn-sincvr_phasn; 
	mask2=(mask>0).*(diffmap<pi).*(diffmap>-pi);
	x_pha(subj,:)=fminsearch(@(x) funfunpha(x,torcvr_phasn(mask2>0),sincvr_phasn(mask2>0)),[1 0],options);

	x_mag(subj,:)=fminsearch(@(x) funfun(x,abs(torcvr_mags(mask>0)),torcvr_magvars(mask>0),abs(sincvr_mags(mask>0)),sincvr_magvars(mask>0)),[1 0],options);

	[corr_mag(:,:,subj) sig_mag(:,:,subj)]=corrcoef(torcvr_mags(mask>0),sincvr_mags(mask>0));
	
	if subj==exsubj
		exsubjtor=abs(torcvr_mags(mask>0));
		exsubjsin=sincvr_mags(mask>0);
		exsubjtorpha=torcvr_phasn(mask2>0);
		exsubjsinpha=sincvr_phasn(mask2>0);
	end
	
	end

	figure;
	plot(exsubjtor.*100,exsubjsin.*100,'.');
	hold on;
	lim=max([get(gca,'xlim'); get(gca,'ylim')]);
	plot([0 2],[0 2],'k-');
	plot([0 2],[0 2]*x_mag(exsubj,1)+x_mag(exsubj,2).*100,'k-');
	ylim([0 2]);
	xlim([0 2]);
	set(gca,'xtick',[0 0.4 0.8 1.2 1.6 2])
	set(gca,'ytick',[0 0.4 0.8 1.2 1.6 2])
	axis square;

	figure;
	bar(x_mag(subj_inc,1));
	hold on;
	plot([0 10],[1 1],'k--');
	axis square
	[h_magps p_magps]=ttest(x_mag(subj_inc,1)-1);
	annotation('textbox',[0.22 0.75 0.5 0.15],'String',['p=' num2str(round(p_magps,2))],'linestyle','none','fontsize',12);
	title('Comparison of slopes between CVR magnitude estimates');
	ylim([0 1.4]);

	figure;
	bar(x_mag(subj_inc,2).*100);
	hold on;
	plot([0 10],[0 0],'k--');
	axis square
	[h_magint p_magint]=ttest(x_mag(subj_inc,2));
	annotation('textbox',[0.22 0.75 0.5 0.15],'String',['p=' num2str(round(p_magint,2))],'linestyle','none','fontsize',12);
	title('Comparison of intercepts between CVR magnitude estimates');

	function f=funfun(x,tor,torvar,sin,sinvar)
		
		f=sum((sin-(x(1).*tor+x(2))).^2./(sinvar+x(1).^2.*torvar));
		
	return;

	function f=funfunpha(x,tor,sin)
			
		f=sum(((x(1).*tor+x(2)-sin)./sqrt(1+x(1).^2)).^2);
		
	return;

	