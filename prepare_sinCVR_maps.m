function bids=prepare_sinCVR_maps(bids,subj)
	
	%using full data
	
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'mean_func']);
	sincvr_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr_pha=atan2(cope1,cope2); %use atan2 for more robustness
	
	bids(subj).func(1).results(1).name='cvr magnitude/phase maps from sinCVR - full data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(1).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-mag'];
	save_avw(sincvr_mag,bids(subj).func(1).results(1).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(1).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-pha'];
	save_avw(sincvr_pha,bids(subj).func(1).results(1).cvr_pha,'f',scales');
	
	%adjust phase values to median phase across brain and unwrap to be comparable with Toronto protocol
	sincvr_phan=sincvr_pha;
	sincvr_phan(sincvr_phan~=0)=sincvr_phan(sincvr_phan~=0)-median(sincvr_pha(sincvr_pha~=0));
	sincvr_phan(sincvr_phan<-pi)=sincvr_phan(sincvr_phan<-pi)+2*pi;
	sincvr_phan(sincvr_phan>pi)=sincvr_phan(sincvr_phan>pi)+2*pi;
	
	save_avw(sincvr_phan,[bids(subj).func(1).results(1).cvr_pha 'n'],'f',scales');
	
	%using 5 mins data
	
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'mean_func']);
	sincvr5min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr5min_pha=atan2(cope1,cope2);	

	bids(subj).func(1).results(2).name='cvr magnitude/phase maps from sinCVR - 5 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(2).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-mag'];
	save_avw(sincvr5min_mag,bids(subj).func(1).results(2).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(2).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-pha'];
	save_avw(sincvr5min_pha,bids(subj).func(1).results(2).cvr_pha,'f',scales');
	
	%adjust phase values to median phase across brain and unwrap to be comparable with Toronto protocol
	sincvr5min_phan=sincvr5min_pha;
	sincvr5min_phan(sincvr5min_phan~=0)=sincvr5min_phan(sincvr5min_phan~=0)-median(sincvr5min_pha(sincvr5min_pha~=0));
	sincvr5min_phan(sincvr5min_phan<-pi)=sincvr5min_phan(sincvr5min_phan<-pi)+2*pi;
	sincvr5min_phan(sincvr5min_phan>pi)=sincvr5min_phan(sincvr5min_phan>pi)+2*pi;
	
	save_avw(sincvr5min_phan,[bids(subj).func(1).results(2).cvr_pha 'n'],'f',scales');
		
	%using 3 mins data
	
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'mean_func']);
	sincvr3min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr3min_pha=atan2(cope1,cope2);
	
	bids(subj).func(1).results(3).name='cvr magnitude/phase maps from sinCVR - 3 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(3).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-mag'];
	save_avw(sincvr3min_mag,bids(subj).func(1).results(3).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(3).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-pha'];
	save_avw(sincvr3min_pha,bids(subj).func(1).results(3).cvr_pha,'f',scales');
	
	%adjust phase values to median phase across brain and unwrap to be comparable with Toronto protocol
	sincvr3min_phan=sincvr3min_pha;
	sincvr3min_phan(sincvr3min_phan~=0)=sincvr3min_phan(sincvr3min_phan~=0)-median(sincvr3min_pha(sincvr3min_pha~=0));
	sincvr3min_phan(sincvr3min_phan<-pi)=sincvr3min_phan(sincvr3min_phan<-pi)+2*pi;
	sincvr3min_phan(sincvr3min_phan>pi)=sincvr3min_phan(sincvr3min_phan>pi)+2*pi;
	
	save_avw(sincvr3min_phan,[bids(subj).func(1).results(3).cvr_pha 'n'],'f',scales');
