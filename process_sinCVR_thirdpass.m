function bids=process_sinCVR_thirdpass(bids,subj)

	%using full data and ev
	
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
	
	bids(subj).func(1).analysis(8).name='petco2 regressors';
	bids(subj).func(1).analysis(8).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev.feat/'];
	bids(subj).func(1).analysis(8).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(8).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',bids(subj).func(1).fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(8).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev.fsf']);
	
	%using 6mins data and ev
	
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_ev6min'];	
	vols=180;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=bold(1:vols);
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
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev6min.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(9).name='petco2 regressors';
	bids(subj).func(1).analysis(9).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev6min.feat/'];
	bids(subj).func(1).analysis(9).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev6min.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev6min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev6min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev6min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(9).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(9).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev6min.fsf']);
	
	%using 5mins data and ev

	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_ev5min'];	
	vols=150;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);	
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=bold(1:vols);
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
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev5min.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(10).name='petco2 regressors';
	bids(subj).func(1).analysis(10).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev5min.feat/'];
	bids(subj).func(1).analysis(10).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev5min.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev5min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev5min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev5min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(10).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(10).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev5min.fsf']);
	
	%using 4mins data and ev
	
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_ev4min'];
	vols=120;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=bold(1:vols);
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
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev4min.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(11).name='petco2 regressors';
	bids(subj).func(1).analysis(11).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev4min.feat/'];
	bids(subj).func(1).analysis(11).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev4min.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev4min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev4min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev4min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(11).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(11).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev4min.fsf']);
	
	%using 3mins data and ev

	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_ev3min'];
	vols=90;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=bold(1:vols);
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
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev3min.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(12).name='petco2 regressors';
	bids(subj).func(1).analysis(12).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev3min.feat/'];
	bids(subj).func(1).analysis(12).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev3min.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev3min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev3min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev3min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(12).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(12).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev3min.fsf']);	
	
	%using 2mins data and ev

	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_ev2min'];
	vols=60;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=bold(1:vols);
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
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev2min.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(13).name='petco2 regressors';
	bids(subj).func(1).analysis(13).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev2min.feat/'];
	bids(subj).func(1).analysis(13).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev2min.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev2min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev2min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev2min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(13).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(13).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev2min.fsf']);

	%using 1mins data and ev
	
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_ev1min'];
	vols=30;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);	
	bold=load([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_braintc.tsv']);
	bold=bold(1:vols);
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
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	[acor lag]=xcorr(petco2i-mean(petco2i),bold-mean(bold));
	[~,I]=max(abs(acor));
	timelag=lag(I)*2;
	ind1=find(bbb.data(:,1)>(events{1}(row)+timelag/60),1,'first');
	ind2=find(bbb.data(:,1)<(events{2}(row)+timelag/60),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1:ind2,2);
	timei=(1:vols)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	petco2ev=petco2i-mean(petco2i);
	petco2ev=petco2ev./max(abs(petco2ev));
	dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev1min.tsv'],petco2ev);
	
	bids(subj).func(1).analysis(14).name='petco2 regressors';
	bids(subj).func(1).analysis(14).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_ev1min.feat/'];
	bids(subj).func(1).analysis(14).ev=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_petco2ev1min.tsv'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev1min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_ev1min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev1min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(14).feat(1:end-1));
		s=strrep(s,'slice_order',bids(subj).func(1).sliceorder);
		s=strrep(s,'func_data',fname);
		s=strrep(s,'fieldmapdata',bids(subj).fmap.fieldmap);
		s=strrep(s,'fieldmapmag',bids(subj).fmap.fieldmap_mag);
		s=strrep(s,'anat_brain',bids(subj).anat.bet);
		s=strrep(s,'petco2_ev',bids(subj).func(1).analysis(14).ev);
		fprintf(fout,'%s\n',s);
		%disp(s)
	end
	fclose(fin);
	fclose(fout);
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_ev1min.fsf']);