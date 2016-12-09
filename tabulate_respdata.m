function T=tabulate_respdata(bids)
	
	dp=1; %decimal places
	
	for subj=1:10
		sincvr_petco2norm(subj,:)=bids(subj).func(1).results(1).petco2base;
		torcvr_petco2norm(subj,:)=bids(subj).func(2).results(1).petco2base;
		sincvr_petco2delta(subj,:)=bids(subj).func(1).results(1).petco2delta;
		torcvr_petco2delta(subj,:)=bids(subj).func(2).results(1).petco2delta;
		subjects{subj,:}=bids(subj).name;
	end
	
	T=table(subjects,round(torcvr_petco2norm,dp),round(torcvr_petco2delta,dp),round(sincvr_petco2norm,dp),round(sincvr_petco2delta,dp));
	
	T.Properties.VariableNames={'Subject' 'TorPetCO2norm' 'TorDeltaPetCO2' 'SinPetCO2norm' 'SinDeltaPetCO2'};
	T.Properties.VariableUnits={'' 'mmHg' 'mmHg' 'mmHg' 'mmHg'};
	
	group(1,:)={'Mean' round(mean(T.TorPetCO2norm),dp) round(mean(T.TorDeltaPetCO2),dp) round(mean(T.SinPetCO2norm),dp) round(mean(T.SinDeltaPetCO2),dp)};
	group(2,:)={'SD' round(std(T.TorPetCO2norm),dp) round(std(T.TorDeltaPetCO2),dp) round(std(T.SinPetCO2norm),dp) round(std(T.SinDeltaPetCO2),dp)};
	
	groupT=cell2table(group);
	groupT.Properties.VariableNames=T.Properties.VariableNames;
	
	T=[T; groupT]
	
	%statistics
	disp('Paired T-test for difference in group mean normocapnic baseline');
	[hnorm pnorm]=ttest(T.TorPetCO2norm(1:10)-T.SinPetCO2norm(1:10)) %are the normocapnic baselines the same
	disp('Paired T-test for difference in group mean change in etCO2');
	[hdelta pdelta]=ttest(T.TorDeltaPetCO2(1:10)-T.SinDeltaPetCO2(1:10)) %are the changes in etCO2 the same
	
	
	