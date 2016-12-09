function T=tabulate_stats(bids);

	dp=2; %decimal places

	for subj=1:10

		%read in data
		[tor_fstat dims scales]=read_avw([bids(subj).func(2).analysis(2).feat 'thresh_zstat1']);
		sin7_fstat=read_avw([bids(subj).func(1).analysis(1).feat 'thresh_zfstat1']);
		sin5_fstat=read_avw([bids(subj).func(1).analysis(3).feat 'thresh_zfstat1']);
		sin3_fstat=read_avw([bids(subj).func(1).analysis(5).feat 'thresh_zfstat1']);
		
		torcvr_mc=load([bids(subj).func(2).analysis(1).feat 'mc/prefiltered_func_data_mcf_rel_mean.rms']);
		sin7cvr_mc=load([bids(subj).func(1).analysis(1).feat 'mc/prefiltered_func_data_mcf_rel_mean.rms']);
		sin5cvr_mc=load([bids(subj).func(1).analysis(3).feat 'mc/prefiltered_func_data_mcf_rel_mean.rms']);
		sin3cvr_mc=load([bids(subj).func(1).analysis(5).feat 'mc/prefiltered_func_data_mcf_rel_mean.rms']);
		
		slicerm=[1 dims(3)];
	
		tor_fstat(:,:,slicerm)=0;
		sin7_fstat(:,:,slicerm)=0;
		sin5_fstat(:,:,slicerm)=0;
		sin3_fstat(:,:,slicerm)=0;
	
		tor_fstats=tor_fstat(:);
		sin7_fstats=sin7_fstat(:);
		sin5_fstats=sin5_fstat(:);
		sin3_fstats=sin3_fstat(:);
	
		tor_n=sum(tor_fstats>0);
		sin7_n=sum(sin7_fstats>0);
		sin5_n=sum(sin5_fstats>0);
		sin3_n=sum(sin3_fstats>0);
	
		n(subj,:)=[tor_n sin7_n sin5_n sin3_n];
		mc(subj,:)=[torcvr_mc sin7cvr_mc sin5cvr_mc sin3cvr_mc];
		subjects{subj,:}=bids(subj).name;
	
	end
	
	mc=round(mc,dp);
	
	T=table(subjects,n(:,1),mc(:,1),n(:,2),mc(:,2),n(:,3),mc(:,3),n(:,4),mc(:,4));
	
	T.Properties.VariableNames={'Subject' 'TorVoxels' 'TorMotion' 'Sin7Voxels' 'Sin7Motion' 'Sin5Voxels' 'Sin5Motion' 'Sin3Voxels' 'Sin3Motion'};
	T.Properties.VariableUnits={'' '' 'mm' '' 'mm' '' 'mm' '' 'mm'};
	
	group(1,:)={'Mean' ceil(mean(T.TorVoxels)) round(mean(T.TorMotion),2) ceil(mean(T.Sin7Voxels)) round(mean(T.Sin7Motion),2) ceil(mean(T.Sin5Voxels)) round(mean(T.Sin5Motion),2) ceil(mean(T.Sin3Voxels)) round(mean(T.Sin3Motion),2)};
	group(2,:)={'SD' ceil(std(T.TorVoxels)) round(std(T.TorMotion),2) ceil(std(T.Sin7Voxels)) round(std(T.Sin7Motion),2) ceil(std(T.Sin5Voxels)) round(std(T.Sin5Motion),2) ceil(std(T.Sin3Voxels)) round(std(T.Sin3Motion),2)};
	
	groupT=cell2table(group);
	groupT.Properties.VariableNames=T.Properties.VariableNames;
	
	T=[T; groupT]
	