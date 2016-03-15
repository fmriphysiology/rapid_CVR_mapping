function bids=prepare_sliceorder(bids,subj)

	for file=1:length(bids(subj).func)
		[timing sliceorder]=sort(bids(subj).func(file).params.SliceTiming');
		%directories created previously
		%if ~exist([bids(subj).dir 'derivatives/' bids(subj).name])
		%	system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name]);
		%end
		%if ~exist([bids(subj).dir 'derivatives/' bids(subj).name '/func'])
		%	system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name '/func']);
		%end
		s=regexp(bids(subj).func(file).fname,'/');
		dlmwrite([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(file).fname(s(end)+1:end) '_sliceorder.tsv'],sliceorder);
		if exist([bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(file).fname(s(end)+1:end) '_sliceorder.tsv'],'file');
			bids(subj).func(file).sliceorder=[bids(subj).dir 'derivatives/' bids(subj).name '/func/' bids(subj).func(file).fname(s(end)+1:end) '_sliceorder.tsv'];
		end
	end