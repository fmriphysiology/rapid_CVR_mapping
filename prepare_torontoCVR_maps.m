function bids=prepare_torontoCVR_maps(bids,subj)

	%from first pass analysis
	
	[cope dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(1).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(1).feat 'mean_func']);
	cvr1_mag=cope./mean_func;
	
	bids(subj).func(2).results(1).name='cvr magnitude map from torontoCVR first pass';
	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(1).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_cvr1-mag'];
	save_avw(cvr1_mag,bids(subj).func(2).results(1).cvr_mag,'f',scales');
	
	%from second pass analysis
	
	[cope dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'mean_func']);
	cvr2_mag=cope./mean_func;
	
	bids(subj).func(2).results(2).name='cvr magnitude map from torontoCVR second pass';
	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(2).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_cvr2-mag'];
	save_avw(cvr2_mag,bids(subj).func(2).results(2).cvr_mag,'f',scales');	
	
	%from second pass - transfer function analysis
	
	ev=dlmread([bids(subj).func(2).analysis(2).feat 'design.mat'],'\t',5,0);
	[func_data func_dims]=read_avw([bids(subj).func(2).analysis(2).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4))';
	
	[cevbold f]=mscohere(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	[pev f]=pwelch(ev(:,1)-mean(ev(:,1)),50,[],[],1/2);
	[pevbold f]=cpsd(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	H=pevbold./repmat(pev,1,64*64*24);
	
	freq=find(f>0.01,1,'first');
	cvr3_mag=reshape(abs(H(freq,:)),func_dims(1),func_dims(2),func_dims(3))./mean_func;
	cvr3_delay=reshape(atan(real(H(freq,:))./imag(H(freq,:))),func_dims(1),func_dims(2),func_dims(3));
	cvr3_delay(cvr3_delay<0)=cvr3_delay(cvr3_delay<0)+pi;

	bids(subj).func(2).results(3).name='cvr magnitude/delay maps from torontoCVR transfer function analysis';
	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(3).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_cvr3-mag'];
	save_avw(cvr3_mag,bids(subj).func(2).results(3).cvr_mag,'f',scales');	
	
	bids(subj).func(2).results(3).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_cvr3-delay'];
	save_avw(cvr3_delay,bids(subj).func(2).results(3).cvr_delay,'f',scales');
	
	