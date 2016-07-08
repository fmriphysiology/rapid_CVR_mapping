function bids=process_torontoCVR_secondpass(bids,subj)

	%analysis using PetCO2 regressor 

	status=system(['fslmeants -i ' bids(subj).func(2).analysis(1).feat 'filtered_func_data -m '...
		bids(subj).func(2).analysis(1).feat 'mask -o '...
		bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_braintc.tsv']);
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_braintc.tsv']);
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
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:210)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:210)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_petco2ev.tsv'],petco2ev);
	
	bids(subj).func(2).analysis(2).name='petco2 regressors';
	bids(subj).func(2).analysis(2).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_bold_ev.feat/'];
	bids(subj).func(2).analysis(2).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_petco2ev.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-torontoCVR_ev.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-torontoCVR_ev.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-torontoCVR_ev.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(2).analysis(2).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(2).sliceorder);
		s=strrep(s,'func_data',bids(subj).func(2).fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(2).analysis(2).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-torontoCVR_ev.fsf']);