function bids=prepare_sinCVR_maps(bids,subj)
	
	%from first pass analysis
	
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'mean_func']);
	cvr1_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	cvr1_delay=atan(cope1./cope2);
	
	bids(subj).func(1).results(1).name='cvr magnitude/delay maps from sinCVR first pass';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(1).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr1-mag'];
	save_avw(cvr1_mag,bids(subj).func(1).results(1).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(1).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr1-delay'];
	save_avw(cvr1_delay,bids(subj).func(1).results(1).cvr_delay,'f',scales');
	
	%from second pass analysis
	%using petco2 regressor
	
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'mean_func']);
	cvr2_mag=cope./mean_func;
	
	bids(subj).func(1).results(2).name='cvr magnitude map from sinCVR - petco2 regressor';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(2).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_cvr2-mag'];
	save_avw(cvr2_mag,bids(subj).func(1).results(2).cvr_mag,'f',scales');
	
	%5 mins data
	
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'mean_func']);
	cvr3_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	cvr3_delay=atan(cope1./cope2);	

	bids(subj).func(1).results(3).name='cvr magnitude/delay maps from sinCVR - 5 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(3).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr3-mag'];
	save_avw(cvr3_mag,bids(subj).func(1).results(3).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(3).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr3-delay'];
	save_avw(cvr3_delay,bids(subj).func(1).results(3).cvr_delay,'f',scales');
		
	%3 mins data
	
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'mean_func']);
	cvr4_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	cvr4_delay=atan(cope1./cope2);
	
	bids(subj).func(1).results(4).name='cvr magnitude/delay maps from sinCVR - 3 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(4).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr4-mag'];
	save_avw(cvr4_mag,bids(subj).func(1).results(4).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(4).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr4-delay'];
	save_avw(cvr4_delay,bids(subj).func(1).results(4).cvr_delay,'f',scales');
	
	%repeat using transfer function analysis
	
	%TFA - 7 mins data
	
	ev=dlmread([bids(subj).func(2).analysis(2).feat 'design.mat'],'\t',5,0);
	[func_data func_dims]=read_avw([bids(subj).func(1).analysis(2).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4))';
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'mean_func']);
	
	[cevbold f]=mscohere(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	[pev f]=pwelch(ev(:,1)-mean(ev(:,1)),50,[],[],1/2);
	[pevbold f]=cpsd(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	H=pevbold./repmat(pev,1,64*64*24);
	
	freq=find(f>(1/60),1,'first');
	cvr5_mag=reshape(abs(H(freq,:)),func_dims(1),func_dims(2),func_dims(3))./mean_func;
	cvr5_delay=reshape(atan(real(H(freq,:))./imag(H(freq,:))),func_dims(1),func_dims(2),func_dims(3));
	cvr5_delay(cvr5_delay<0)=cvr5_delay(cvr5_delay<0)+pi;
	
	bids(subj).func(1).results(5).name='cvr magnitude/delay maps from sinCVR - transfer function analysis';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(5).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr5-mag'];
	save_avw(cvr5_mag,bids(subj).func(1).results(5).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(5).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr5-delay'];
	save_avw(cvr5_delay,bids(subj).func(1).results(5).cvr_delay,'f',scales');
	
	%TFA - 5 mins data
	
	ev=ev(1:5*60/2,:);
	[func_data func_dims]=read_avw([bids(subj).func(1).analysis(3).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4))';
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'mean_func']);
	
	[cevbold f]=mscohere(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	[pev f]=pwelch(ev(:,1)-mean(ev(:,1)),50,[],[],1/2);
	[pevbold f]=cpsd(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	H=pevbold./repmat(pev,1,64*64*24);
	
	freq=find(f>(1/60),1,'first');
	cvr6_mag=reshape(abs(H(freq,:)),func_dims(1),func_dims(2),func_dims(3))./mean_func;
	cvr6_delay=reshape(atan(real(H(freq,:))./imag(H(freq,:))),func_dims(1),func_dims(2),func_dims(3));
	cvr6_delay(cvr6_delay<0)=cvr6_delay(cvr6_delay<0)+pi;
	
	bids(subj).func(1).results(6).name='cvr magnitude/delay maps from sinCVR - transfer function analysis 5 mins';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(6).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr6-mag'];
	save_avw(cvr6_mag,bids(subj).func(1).results(6).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(6).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr6-delay'];
	save_avw(cvr6_delay,bids(subj).func(1).results(6).cvr_delay,'f',scales');

	%TFA - 3 mins data
	
	ev=ev(1:3*60/2,:);
	[func_data func_dims]=read_avw([bids(subj).func(1).analysis(4).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4))';
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'mean_func']);
	
	[cevbold f]=mscohere(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	[pev f]=pwelch(ev(:,1)-mean(ev(:,1)),50,[],[],1/2);
	[pevbold f]=cpsd(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	H=pevbold./repmat(pev,1,64*64*24);
	
	freq=find(f>(1/60),1,'first');
	cvr7_mag=reshape(abs(H(freq,:)),func_dims(1),func_dims(2),func_dims(3))./mean_func;
	cvr7_delay=reshape(atan(real(H(freq,:))./imag(H(freq,:))),func_dims(1),func_dims(2),func_dims(3));
	cvr7_delay(cvr7_delay<0)=cvr7_delay(cvr7_delay<0)+pi;
	
	bids(subj).func(1).results(7).name='cvr magnitude/delay maps from sinCVR - transfer function analysis 3 mins';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(7).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr7-mag'];
	save_avw(cvr7_mag,bids(subj).func(1).results(7).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(7).cvr_delay=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_cvr7-delay'];
	save_avw(cvr7_delay,bids(subj).func(1).results(7).cvr_delay,'f',scales');
	
	