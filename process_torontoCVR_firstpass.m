function bids=process_torontoCVR_firstpass(bids,subj)

	%full analysis using boxcar model rather than PetCO2 regressor

	bids(subj).func(2).analysis(1).name='boxcar regressors';
	bids(subj).func(2).analysis(1).feat=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).name '_task-torontoCVR_bold.feat/'];
	system(['cp ' bids(subj).dir 'derivatives/code/feat_designs/task-torontoCVR.* ' bids(subj).dir 'derivatives/' bids(subj).name '/func/']);
	fin=fopen([bids(subj).dir 'derivatives/code/feat_designs/task-torontoCVR.fsf'],'r');
	fout=fopen([bids(subj).dir 'derivatives/' bids(subj).name '/func/task-torontoCVR.fsf'],'w');
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
	status=system(['feat ' bids(subj).dir 'derivatives/' bids(subj).name '/func/task-torontoCVR.fsf']);