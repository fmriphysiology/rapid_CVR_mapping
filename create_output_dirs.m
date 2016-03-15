function bids=create_output_dirs(bids,subj);

		%create subj dir
		if ~exist([bids(subj).dir 'derivatives/' bids(subj).name])
			system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name]);
		end
		
		%create subj/func dir
		if ~exist([bids(subj).dir 'derivatives/' bids(subj).name '/func'])
			system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name '/func']);
		end
		
		%create subj/fmap dir
		if ~exist([bids(subj).dir 'derivatives/' bids(subj).name '/fmap'])
			system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name '/fmap']);
		end
		
		%create subj/anat dir
		if ~exist([bids(subj).dir 'derivatives/' bids(subj).name '/anat'])
			system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name '/anat']);
		end
		
		%create subj/results dir
		if ~exist([bids(subj).dir 'derivatives/' bids(subj).name '/results'])
			system(['mkdir ' bids(subj).dir 'derivatives/' bids(subj).name '/results']);
		end