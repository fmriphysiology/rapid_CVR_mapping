function bids=process_sinCVR_secondpass(bids,subj)

	%reducing the amount of data used in the analysis 
	
	% 6 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_6min'];
	vols=180;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(2).name='sin/cos regressors - 6 mins data';
	bids(subj).func(1).analysis(2).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_6mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_6min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_6min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_6min.fsf'],'w');
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
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_6min.fsf']);
	
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

	% 4 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_4min'];
	vols=120;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(4).name='sin/cos regressors - 4 mins data';
	bids(subj).func(1).analysis(4).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_4mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_4min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_4min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_4min.fsf'],'w');
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
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_4min.fsf']);

	% 3 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_3min'];
	vols=90;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(5).name='sin/cos regressors - 3 mins data';
	bids(subj).func(1).analysis(5).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_3mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_3min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_3min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_3min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(5).feat(1:end-1));
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
	
	% 2 mins of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_2min'];
	vols=60;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(6).name='sin/cos regressors - 2 mins data';
	bids(subj).func(1).analysis(6).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_2mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_2min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_2min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_2min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(6).feat(1:end-1));
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
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_2min.fsf']);
	
	% 1 min of data
	s=regexp(bids(subj).func(1).fname,'/');
	fname=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(1).fname(s(end)+1:end) '_1min'];
	vols=30;
	status=system(['fslroi ' bids(subj).func(1).fname ' ' fname ' 0 ' num2str(vols)]);
	bids(subj).func(1).analysis(7).name='sin/cos regressors - 1 mins data';
	bids(subj).func(1).analysis(7).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-sinCVR_bold_1mins.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_1min.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-sinCVR_1min.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_1min.fsf'],'w');
	while ~feof(fin)
		s=fgetl(fin);
		s=strrep(s,'output_dir',bids(subj).func(1).analysis(7).feat(1:end-1));
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
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-sinCVR_1min.fsf']);