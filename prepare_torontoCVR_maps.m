function bids=prepare_torontoCVR_maps(bids,subj)
	
	%cvr magnitude from regression analysis 
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(2).analysis(2).feat 'stats/sigmasquareds']);
	
	%read in resp data and estimate change in PetCO2
	events_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_events.tsv'];
	fid=fopen(events_in,'r');
	events=textscan(fid,'%f %f %s','delimiter','\t','headerlines',1);
	fclose(fid);
	
	row=find(strcmp(events{3},'torontoCVR'));
	bbb_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_bbb.tsv'];
	bbb=importdata(bbb_in,'\t');
	ind1=find(bbb.data(:,1)>events{1}(row),1,'first');
	ind2=find(bbb.data(:,1)<events{2}(row),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	torpetco2=bbb.data(ind1:ind2,2);
	timei=(1:210)'.*2;
	torstim=(timei>60).*(timei<(60+45))+(timei>(60+45+90)).*(timei<(60+45+90+120));
	torpetco2i=interp1(time,torpetco2,timei,'linear','extrap');
	[acor lag]=xcorr(torpetco2i-mean(torpetco2i),torstim-mean(torstim));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	torpetco2=bbb.data(ind1:ind2,2);
	timei=(1:210)'.*2;
	torpetco2i=interp1(time,torpetco2,timei,'linear','extrap');
	torX=[torstim ones(size(torstim))];
	[tora torstd]=lscov(torX,torpetco2i);
	tordeltapetco2=tora(1);
	torvardeltapetco2=torstd(1)^2;
			
	bids(subj).func(2).results(1).petco2base=tora(2);
	bids(subj).func(2).results(1).petco2basevar=torstd(2);
	bids(subj).func(2).results(1).petco2delta=tordeltapetco2;
	bids(subj).func(2).results(1).petco2deltavar=torvardeltapetco2;

	%estimate magnitude and errors	
	torcvr_mag=cope./mean_func./tordeltapetco2;
	torcvr_mag(isnan(torcvr_mag))=0;
	
	torcvr_magvar=(cope./mean_func./tordeltapetco2).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2+torvardeltapetco2./tordeltapetco2.^2);
		
	%write out mag images
	bids(subj).func(2).results(1).name='cvr magnitude/delay map from torontoCVR';
	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(1).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_torcvr-mag'];
	save_avw(torcvr_mag,bids(subj).func(2).results(1).cvr_mag,'f',scales');	
	
	bids(subj).func(2).results(1).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_torcvr-magvar'];
	save_avw(torcvr_magvar,bids(subj).func(2).results(1).cvr_magvar,'f',scales');	
	
	%cvr phase from transfer function analysis (no errors since no defined method for calculating variance)
	ev=dlmread([bids(subj).func(2).analysis(2).feat 'design.mat'],'\t',5,0);
	[func_data func_dims]=read_avw([bids(subj).func(2).analysis(2).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4))';
	
	[H F]=tfestimate(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	freq=find(F>0.01,1,'first');
	torcvr_pha=reshape(atan2(imag(H(freq,:)),real(H(freq,:))),func_dims(1),func_dims(2),func_dims(3)); %use atan2 for more robustness

	s=regexp(bids(subj).func(2).fname,'/');
	bids(subj).func(2).results(1).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(2).fname(s(end)+1:end) '_torcvr-pha'];
	save_avw(torcvr_pha,bids(subj).func(2).results(1).cvr_pha,'f',scales');
	
	