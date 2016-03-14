for subj=9:10;

	bids_dir='/Users/nickb/Analysis/fmrib/cvr_study/';

	% initialise bids data structure
	bids=bids_init(bids_dir);

	% prepare fieldmaps
	deltaTE=(str2num(bids(subj).fmap.params.EchoTime2)-str2num(bids(subj).fmap.params.EchoTime1))*1000;
	status=system([bids_dir 'derivatives/code/prep_fmap.sh ' bids_dir ' ' bids(subj).name ' ' num2str(deltaTE)]);
	if status==0
		bids(subj).fmap.fieldmap=[bids_dir 'derivatives/' bids(subj).name '/fmap/' bids(subj).name '_fieldmap'];
		bids(subj).fmap.fieldmap_mag=[bids_dir 'derivatives/' bids(subj).name '/fmap/' bids(subj).name '_magnitude_e2_brain_ero'];
	end

	% prepare anatomical images
	status=system([bids_dir 'derivatives/code/prep_anat.sh ' bids_dir ' ' bids(subj).name]);
	if status==0
		bids(subj).anat.bet=[bids_dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain'];
	end

	% prepare slice order file
	for file=1:length(bids(subj).func)
		[timing sliceorder]=sort(bids(subj).func(file).params.SliceTiming');
		if ~exist([bids_dir 'derivatives/' bids(subj).name])
			system(['mkdir ' bids_dir 'derivatives/' bids(subj).name]);
		end
		if ~exist([bids_dir 'derivatives/' bids(subj).name '/func'])
			system(['mkdir ' bids_dir 'derivatives/' bids(subj).name '/func']);
		end
		s=regexp(bids(subj).func(file).fname,'/');
		dlmwrite([bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(file).fname(s(end)+1:end) '_sliceorder.tsv'],sliceorder);
		if exist([bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(file).fname(s(end)+1:end) '_sliceorder.tsv'],'file');
			bids(subj).func(file).sliceorder=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(file).fname(s(end)+1:end) '_sliceorder.tsv'];
		end
	end

	% prepare 1st pass sinCVR feat design
	bids(subj).func(1).analysis(1).name='sin/cos regressors';
	bids(subj).func(1).analysis(1).feat=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold.feat/'];
	system(['cp ' bids_dir 'derivatives/code/feat_designs/task-sinCVR.* ' bids_dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids_dir 'derivatives/code/feat_designs/task-sinCVR.fsf'],'r');
	fout=fopen([bids_dir 'derivatives/' bids(subj).name '/func/task-sinCVR.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(1).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',bids(subj).func(1).fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids_dir 'derivatives/' bids(subj).name '/func/task-sinCVR.fsf']);

	% prepare 1st pass torontoCVR feat design
	bids(subj).func(2).analysis(1).name='boxcar regressors';
	bids(subj).func(2).analysis(1).feat=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_bold.feat/'];
	system(['cp ' bids_dir 'derivatives/code/feat_designs/task-torontoCVR.* ' bids_dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids_dir 'derivatives/code/feat_designs/task-torontoCVR.fsf'],'r');
	fout=fopen([bids_dir 'derivatives/' bids(subj).name '/func/task-torontoCVR.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(2).analysis(1).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(2).sliceorder);
		s=strrep(s,'func_data',bids(subj).func(2).fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids_dir 'derivatives/' bids(subj).name '/func/task-torontoCVR.fsf']);
	
	% prepare 2nd pass sinCVR feat designs - data reduction
	% 5 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_5min'];
	vols=150;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(2).name='sin/cos regressors - 5 mins data';
	bids(subj).func(1).analysis(2).feat=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_5mins.feat/'];
	system(['cp ' bids_dir 'derivatives/code/feat_designs/task-sinCVR_5min.* ' bids_dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids_dir 'derivatives/code/feat_designs/task-sinCVR_5min.fsf'],'r');
	fout=fopen([bids_dir 'derivatives/' bids(subj).name '/func/task-sinCVR_5min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(2).feat(1:end-1));
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
	status=system(['feat ' bids_dir 'derivatives/' bids(subj).name '/func/task-sinCVR_5min.fsf']);

	% 3 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_3min'];
	vols=90;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(3).name='sin/cos regressors - 3 mins data';
	bids(subj).func(1).analysis(3).feat=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_3mins.feat/'];
	system(['cp ' bids_dir 'derivatives/code/feat_designs/task-sinCVR_3min.* ' bids_dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids_dir 'derivatives/code/feat_designs/task-sinCVR_3min.fsf'],'r');
	fout=fopen([bids_dir 'derivatives/' bids(subj).name '/func/task-sinCVR_3min.fsf'],'w');
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
	status=system(['feat ' bids_dir 'derivatives/' bids(subj).name '/func/task-sinCVR_3min.fsf']);

	% prepare 2nd pass torontoCVR feat designs - subject specific PetCO2 regressor
	%status=system(['fslmaths ' bids(subj).func(2).analysis(1).feat 'mean_func -bin '...
		bids(subj).func(2).analysis(1).feat 'bet_mask']); %mask already exists
	status=system(['fslmeants -i ' bids(subj).func(2).analysis(1).feat 'filtered_func_data -m '...
		bids(subj).func(2).analysis(1).feat 'mask -o '...
		bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_braintc.tsv']);
	bold=load([bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_braintc.tsv']);
	events_in=[bids_dir bids(subj).name '/resp/' bids(subj).name '_respdata_events.tsv'];
	fid=fopen(events_in,'r');
	events=textscan(fid,'%f %f %s','delimiter','\t','headerlines',1);
	fclose(fid);
	row=find(strcmp(events{3},'torontoCVR'));
	bbb_in=[bids_dir bids(subj).name '/resp/' bids(subj).name '_respdata_bbb.tsv'];
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
	dlmwrite([bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_petco2ev.tsv'],petco2ev);
	
	bids(subj).func(2).analysis(2).name='petco2 regressors';
	bids(subj).func(2).analysis(2).feat=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_bold_ev.feat/'];
	bids(subj).func(2).analysis(2).ev=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_petco2ev.tsv'];
	system(['cp ' bids_dir 'derivatives/code/feat_designs/task-torontoCVR_ev.* ' bids_dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids_dir 'derivatives/code/feat_designs/task-torontoCVR_ev.fsf'],'r');
	fout=fopen([bids_dir 'derivatives/' bids(subj).name '/func/task-torontoCVR_ev.fsf'],'w');
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
	status=system(['feat ' bids_dir 'derivatives/' bids(subj).name '/func/task-torontoCVR_ev.fsf']);



end