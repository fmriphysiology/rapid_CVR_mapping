function figure_compare_tor_sin(bids)

	for subj=1:10

	%read in toronto CVR maps 
	[torcvr_mag dims scales]=read_avw(bids(subj).func(2).results(1).cvr_mag);
	torcvr_magvar=read_avw(bids(subj).func(2).results(1).cvr_magvar);
	torcvr_pha=read_avw(bids(subj).func(2).results(1).cvr_pha);

	%read in sinusoid CVR maps (full data)
	sincvr_mag=read_avw([bids(subj).func(1).results(1).cvr_mag '_reg']);
	sincvr_magvar=read_avw([bids(subj).func(1).results(1).cvr_magvar '_reg']);
	sincvr_pha=read_avw([bids(subj).func(1).results(1).cvr_pha '_reg']);
	sincvr_phavar=read_avw([bids(subj).func(1).results(1).cvr_phavar '_reg']);
	
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
	
	torcvr_mag(:,:,slicerm)=0;
	torcvr_magvar(:,:,slicerm)=0;
	torcvr_pha(:,:,slicerm)=0;
	
	sincvr_mag(:,:,slicerm)=0;
	sincvr_magvar(:,:,slicerm)=0;
	sincvr_pha(:,:,slicerm)=0;
	sincvr_phavar(:,:,slicerm)=0;
	
	%rearrange matrices
	gms=gm(:);
	wms=wm(:);
	atlass=atlas(:);
	
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

	torcvr_magavg(subj,:)=mean(abs(torcvr_mags(mask>0)));
	torcvr_magavg_sd(subj,:)=std(abs(torcvr_mags(mask>0)));
	torcvr_magavg_se(subj,:)=std(abs(torcvr_mags(mask>0)))./sqrt(sum(mask>0));
	sincvr_magavg(subj,:)=mean(sincvr_mags(mask>0));
	sincvr_magavg_sd(subj,:)=std(sincvr_mags(mask>0));
	sincvr_magavg_se(subj,:)=std(sincvr_mags(mask>0))./sqrt(sum(mask>0));
	
	torcvr_phaavg(subj,:)=mean(abs(torcvr_phasn(mask>0)));
	torcvr_phaavg_sd(subj,:)=std(abs(torcvr_phasn(mask>0)));
	torcvr_phaavg_se(subj,:)=std(abs(torcvr_phasn(mask>0)))./sqrt(sum(mask>0));
	sincvr_phaavg(subj,:)=mean(sincvr_phasn(mask>0));
	sincvr_phaavg_sd(subj,:)=std(sincvr_phasn(mask>0));
	sincvr_phaavg_se(subj,:)=std(sincvr_phasn(mask>0))./sqrt(sum(mask>0));
	
	for k=1:9
		torcvr_magroi(k,subj)=mean(abs(100.*torcvr_mags(find((mask>0).*(atlass==k)))));
		torcvr_magroisd(k,subj)=std(abs(100.*torcvr_mags(find((mask>0).*(atlass==k)))));
		torcvr_magroise(k,subj)=std(abs(100.*torcvr_mags(find((mask>0).*(atlass==k)))))./sqrt(sum((mask>0).*(atlass==k)));
		sincvr_magroi(k,subj)=mean(100.*sincvr_mags(find((mask>0).*(atlass==k))));
		sincvr_magroisd(k,subj)=std(100.*sincvr_mags(find((mask>0).*(atlass==k))));
		sincvr_magroise(k,subj)=std(100.*sincvr_mags(find((mask>0).*(atlass==k))))./sqrt(sum((mask>0).*(atlass==k)));
	end
	
	for k=1:9
		torcvr_pharoi(k,subj)=mean(torcvr_phasn(find((mask>0).*(atlass==k))));
		torcvr_pharoisd(k,subj)=std(torcvr_phasn(find((mask>0).*(atlass==k))));
		torcvr_pharoise(k,subj)=std(torcvr_phasn(find((mask>0).*(atlass==k))))./sqrt(sum((mask>0).*(atlass==k)));
		sincvr_pharoi(k,subj)=mean(sincvr_phasn(find((mask>0).*(atlass==k))));
		sincvr_pharoisd(k,subj)=std(sincvr_phasn(find((mask>0).*(atlass==k))));
		sincvr_pharoise(k,subj)=std(sincvr_phasn(find((mask>0).*(atlass==k))))./sqrt(sum((mask>0).*(atlass==k)));
	end
	
	X=[abs(torcvr_mags(mask>0)) ones(size(abs(torcvr_mags(mask>0))))];
	a_mag(:,subj)=X\sincvr_mags(mask>0);
	
	end

	%subjects to include in analysis
	subj_inc=1:9; 

	addpath /Users/nickb/Documents/MATLAB/errorbarxy/
	figure;
	clf;
	nums=[2 3 4 5 6 8];
	c='rgbcmy';
	for k=1:6
		errorbarxy(torcvr_pharoi(nums(k),subj_inc),sincvr_pharoi(nums(k),subj_inc),torcvr_pharoise(nums(k),subj_inc),sincvr_pharoise(nums(k),subj_inc),{[c(k) '.'] c(k) c(k)});
		hold on;
	end

	[r_pha p_pha]=corrcoef(reshape(torcvr_pharoi(nums,subj_inc),[],1),reshape(sincvr_pharoi(nums,subj_inc),[],1));
	ptext='';
	if p_pha(1,2)<0.001
		ptext='p<0.001, ';
	elseif p_pha(1,2)<0.01 
		ptext='p<0.01, ';
	elseif p_pha(1,2)<0.05 
		ptext='p<0.05, ';
	else 
		ptext='';
	end
	[agm_pha agm_pha_sd]=lscov([reshape(torcvr_pharoi(nums,subj_inc),[],1) ones(length(nums)*length(subj_inc),1)],reshape(sincvr_pharoi(nums,subj_inc),[],1));
	axis square;
	lims(1)=min([xlim ylim]);
	lims(2)=max([xlim ylim]);
	xlim(lims);
	ylim(lims);
	legend('Cerebellum','Frontal','Insula','Occipital','Parietal','Temporal','location','southeast')
	fit_pha=fit(reshape(torcvr_pharoi(nums,subj_inc),[],1),reshape(sincvr_pharoi(nums,subj_inc),[],1),'poly1');
	fit_pha_ci_vals=confint(fit_pha,0.95);
	fit_pha_ci=predint(fit_pha,linspace(lims(1),lims(2),100),0.95,'functional','off');
	plot(lims',[lims' ones(size(lims'))]*agm_pha,'k');
	plot(linspace(lims(1),lims(2),100),fit_pha_ci,'k--');
	annotation('textbox',[0.22 0.75 0.5 0.15],'String',['R=' num2str(round(r_pha(1,2),2)) ', ' ptext 'Slope=' num2str(round(agm_pha(1),2)) '+/-' num2str(round(agm_pha_sd(1),2))],'linestyle','none','fontsize',12);
	title('Comparison of CVR phase estimates');
		
	figure;
	clf;
	nums=[2 3 4 5 6 8];
	c='rgbcmy';
	for k=1:6
		errorbarxy(torcvr_magroi(nums(k),subj_inc),sincvr_magroi(nums(k),subj_inc),torcvr_magroise(nums(k),subj_inc),sincvr_magroise(nums(k),subj_inc),{[c(k) '.'] c(k) c(k)});
		hold on;
	end
	
	[r_mag p_mag]=corrcoef(reshape(torcvr_magroi(nums,subj_inc),[],1),reshape(sincvr_magroi(nums,subj_inc),[],1));
	ptext='';
	if p_mag(1,2)<0.001
		ptext='p<0.001, ';
	elseif p_mag(1,2)<0.01 
		ptext='p<0.01, ';
	elseif p_mag(1,2)<0.05 
		ptext='p<0.05, ';
	else 
		ptext='';
	end
	[agm_mag agm_mag_sd]=lscov([reshape(torcvr_magroi(nums,subj_inc),[],1) ones(length(nums)*length(subj_inc),1)],reshape(sincvr_magroi(nums,subj_inc),[],1));
	axis square;
	lims(1)=0;
	lims(2)=max([xlim ylim]);
	xlim(lims);
	ylim(lims);
	legend('Cerebellum','Frontal','Insula','Occipital','Parietal','Temporal','location','southeast')
	fit_mag=fit(reshape(torcvr_magroi(nums,subj_inc),[],1),reshape(sincvr_magroi(nums,subj_inc),[],1),'poly1');
	fit_mag_ci_vals=confint(fit_mag,0.95);
	fit_mag_ci=predint(fit_mag,linspace(lims(1),lims(2),100),0.95,'functional','off');
	plot(lims',[lims' ones(size(lims'))]*agm_mag,'k');
	plot(linspace(lims(1),lims(2),100),fit_mag_ci,'k--');
	annotation('textbox',[0.22 0.75 0.5 0.15],'String',['R=' num2str(round(r_mag(1,2),2)) ', ' ptext 'Slope=' num2str(round(agm_mag(1),2)) '+/-' num2str(round(agm_mag_sd(1),2))],'linestyle','none','fontsize',12);
	title('Comparison of normalised CVR magnitude estimates');
	
	for subj=1:10
		[agm_magps(:,subj) agm_mag_sdps(:,subj)]=lscov([torcvr_magroi(nums,subj) ones(6,1)],sincvr_magroi(nums,subj));
		[agm_phaps(:,subj) agm_pha_sdps(:,subj)]=lscov([torcvr_pharoi(nums,subj) ones(6,1)],sincvr_pharoi(nums,subj));
	end
	
	figure;
	bar(agm_magps(1,subj_inc));
	hold on;
	plot([0 10],[1 1],'k--');
	axis square
	[h_magps p_magps]=ttest(agm_magps(1,subj_inc)-1)
	annotation('textbox',[0.22 0.75 0.5 0.15],'String',['p=' num2str(round(p_magps,2))],'linestyle','none','fontsize',12);
	title('Comparison of slopes between CVR magnitude estimates');

	
	figure;
	bar(agm_phaps(1,subj_inc));
	hold on;
	plot([0 10],[1 1],'k--');	
	axis square
	[h_phaps p_phaps]=ttest(agm_phaps(1,subj_inc)-1)
	annotation('textbox',[0.22 0.75 0.5 0.15],'String',['p=' num2str(round(p_phaps,2))],'linestyle','none','fontsize',12);
	title('Comparison of slopes between CVR phase estimates');
		