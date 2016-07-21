function bids=compare_cvr_maps(bids,subj);

	%read in toronto CVR maps 
	[torcvr_mag dims scales]=read_avw(bids(subj).func(2).results(1).cvr_mag);
	torcvr_magvar=read_avw(bids(subj).func(2).results(1).cvr_magvar);
	torcvr_pha=read_avw(bids(subj).func(2).results(1).cvr_pha);

	%read in sinusoid CVR maps (full data)
	sincvr_mag=read_avw([bids(subj).func(1).results(1).cvr_mag '_reg']);
	sincvr_magvar=read_avw([bids(subj).func(1).results(1).cvr_magvar '_reg']);
	sincvr_pha=read_avw([bids(subj).func(1).results(1).cvr_pha '_reg']);
	sincvr_phavar=read_avw([bids(subj).func(1).results(1).cvr_phavar '_reg']);

	%read in sinusoid CVR maps (5 mins data)
	sincvr5min_mag=read_avw([bids(subj).func(1).results(2).cvr_mag '_reg']);
	sincvr5min_magvar=read_avw([bids(subj).func(1).results(2).cvr_magvar '_reg']);
	sincvr5min_pha=read_avw([bids(subj).func(1).results(2).cvr_pha '_reg']);
	sincvr5min_phavar=read_avw([bids(subj).func(1).results(2).cvr_phavar '_reg']);
	
	%read in sinusoid CVR maps (3 mins data)
	sincvr3min_mag=read_avw([bids(subj).func(1).results(3).cvr_mag '_reg']);
	sincvr3min_magvar=read_avw([bids(subj).func(1).results(3).cvr_magvar '_reg']);
	sincvr3min_pha=read_avw([bids(subj).func(1).results(3).cvr_pha '_reg']);
	sincvr3min_phavar=read_avw([bids(subj).func(1).results(3).cvr_phavar '_reg']);
	
	%read in tissue masks
	gm=read_avw([bids(subj).anat.gm '_reg']);
	wm=read_avw([bids(subj).anat.wm '_reg']);
	
	%discard top and bottom slices (lost to realignment)
	slicerm=[1 dims(3)];
		
	gm(:,:,slicerm)=0;
	wm(:,:,slicerm)=0;
	
	torcvr_mag(:,:,slicerm)=0;
	torcvr_magvar(:,:,slicerm)=0;
	torcvr_pha(:,:,slicerm)=0;
	
	sincvr_mag(:,:,slicerm)=0;
	sincvr_magvar(:,:,slicerm)=0;
	sincvr_pha(:,:,slicerm)=0;
	sincvr_phavar(:,:,slicerm)=0;
	
	sincvr5min_mag(:,:,slicerm)=0;
	sincvr5min_magvar(:,:,slicerm)=0;
	sincvr5min_pha(:,:,slicerm)=0;
	sincvr5min_phavar(:,:,slicerm)=0;
	
	sincvr3min_mag(:,:,slicerm)=0;
	sincvr3min_magvar(:,:,slicerm)=0;
	sincvr3min_pha(:,:,slicerm)=0;
	sincvr3min_phavar(:,:,slicerm)=0;
	
	%rearrange matrices
	gms=gm(:);
	wms=wm(:);
	
	torcvr_mags=torcvr_mag(:);
	torcvr_magvars=torcvr_magvar(:);
	torcvr_phas=torcvr_pha(:);
	
	sincvr_mags=sincvr_mag(:);
	sincvr_magvars=sincvr_magvar(:);
	sincvr_phas=sincvr_pha(:);
	sincvr_phavars=sincvr_phavar(:);
	
	sincvr5min_mags=sincvr5min_mag(:);
	sincvr5min_magvars=sincvr5min_magvar(:);
	sincvr5min_phas=sincvr5min_pha(:);
	sincvr5min_phavars=sincvr5min_phavar(:);
	
	sincvr3min_mags=sincvr3min_mag(:);
	sincvr3min_magvars=sincvr3min_magvar(:);
	sincvr3min_phas=sincvr3min_pha(:);
	sincvr3min_phavars=sincvr3min_phavar(:);
	
	%compare cvr magnitude between methods
	mask=(sincvr_mags~=0).*(torcvr_mags~=0);
	
	options=optimset('fminsearch');
	options.Display='iter';
	
	%toronto vs sinusoid
	x_torsin=fminsearch(@(x) modchisqd(x,sincvr_mags(mask>0),sincvr_magvars(mask>0),abs(torcvr_mags(mask>0)),torcvr_magvars(mask>0)),[1 0]);
	
	%sinusoid vs sinusoid 5min
	x_sin5minsin=fminsearch(@(x) modchisqd(x,sincvr_mags(mask>0),sincvr_magvars(mask>0),sincvr5min_mags(mask>0),sincvr5min_magvars(mask>0)),[1 0]);
	
	%sinusoid vs sinusoid 3min
	x_sin3minsin=fminsearch(@(x) modchisqd(x,sincvr_mags(mask>0),sincvr_magvars(mask>0),sincvr3min_mags(mask>0),sincvr3min_magvars(mask>0)),[1 0]);
	
	%compare cvr phase with sinusoid data reduction
	xpha_sin5minsin=fminsearch(@(x) modchisqd(x,sincvr_phas(mask>0),sincvr_phavars(mask>0),sincvr5min_phas(mask>0),sincvr5min_phavars(mask>0)),[1 0]);
	xpha_sin3minsin=fminsearch(@(x) modchisqd(x,sincvr_phas(mask>0),sincvr_phavars(mask>0),sincvr3min_phas(mask>0),sincvr3min_phavars(mask>0)),[1 0]);

	%save out values
	bids(subj).comparison.tor_vs_sin_mag=x_torsin;
	bids(subj).comparison.sin5min_vs_sin_mag=x_sin5minsin;
	bids(subj).comparison.sin3min_vs_sin_mag=x_sin3minsin;
	bids(subj).comparison.sin5min_vs_sin_pha=xpha_sin5minsin;
	bids(subj).comparison.sin3min_vs_sin_pha=xpha_sin3minsin;	
	
	%keyboard;
	
	
	function f=modchisqd(x,sinvals,sinsigma,torvals,torsigma)
	
		f=sum((torvals-(x(1).*sinvals+x(2))).^2./(torsigma+x(1)^2.*sinsigma));
	
	return;
