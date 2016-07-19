function bids=prepare_torontoCVR_maps(bids,subj)
	
	%cvr magnitude from regression analysis 
	
	[cope dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'mean_func']);
	torcvr_mag=cope./mean_func;
	torcvr_mag(isnan(torcvr_mag))=0;
	
	bids(subj).func(2).results(1).name='cvr magnitude/delay map from torontoCVR';
	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(1).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_torcvr-mag'];
	save_avw(torcvr_mag,bids(subj).func(2).results(1).cvr_mag,'f',scales');	
	
	%cvr phase from transfer function analysis
	
	ev=dlmread([bids(subj).func(2).analysis(2).feat 'design.mat'],'\t',5,0);
	[func_data func_dims]=read_avw([bids(subj).func(2).analysis(2).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4))';
	
	[H F]=tfestimate(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	freq=find(F>0.01,1,'first');
	torcvr_pha=reshape(atan2(imag(H(freq,:)),real(H(freq,:))),func_dims(1),func_dims(2),func_dims(3)); %use atan2 for more robustness

	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(1).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_torcvr-pha'];
	save_avw(torcvr_pha,bids(subj).func(2).results(1).cvr_pha,'f',scales');
	
	