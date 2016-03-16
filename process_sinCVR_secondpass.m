function bids=process_sinCVR_secondpass(bids,subj)

	% using petco2 regressor
	status=system(['fslmeants -i ' bids(subj).func(1).analysis(1).feat 'filtered_func_data -m '...
		bids(subj).func(1).analysis(1).feat 'mask -o '...
		bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	events_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_events.tsv'];
	fid=fopen(events_in,'r');
	events=textscan(fid,'%f %f %s','delimiter','\t','headerlines',1);
	fclose(fid);
	row=find(strcmp(events{3},'sinCVR'));
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
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(2).name='petco2 regressors';
	bids(subj).func(1).analysis(2).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev.feat/'];
	bids(subj).func(1).analysis(2).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(2).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',bids(subj).func(1).fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(2).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev.fsf']);

	% 5 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_5min'];
	vols=150;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(3).name='sin/cos regressors - 5 mins data';
	bids(subj).func(1).analysis(3).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_5mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_5min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_5min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_5min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(3).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_5min.fsf']);

	% 3 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_3min'];
	vols=90;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(4).name='sin/cos regressors - 3 mins data';
	bids(subj).func(1).analysis(4).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_3mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_3min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_3min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_3min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(4).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_3min.fsf']);
	

	
	