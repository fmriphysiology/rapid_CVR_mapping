function bids=bids_init(bids_dir)

if bids_dir(end)~='/'
	bids_dir=[bids_dir '/'];
end

participants=importdata([bids_dir 'participants.tsv']);

num_subj=size(participants.data,1);

anat(1).fname='_T1w';
fmap(1).magnitude='_magnitude';
fmap(1).phasediff='_phasediff';
func(1).fname='_task-sinCVR_bold';
func(2).fname='_task-torontoCVR_bold';

% setup data structure
for subj=1:num_subj

	bids(subj).name=participants.textdata{subj+1};
	
	for file=1:length(anat)
		bids(subj).anat(file).fname=[bids_dir bids(subj).name '/anat/' bids(subj).name anat(file).fname];
	end
	
	for file=1:length(fmap)
		bids(subj).fmap(file).magnitude=[bids_dir bids(subj).name '/fmap/' bids(subj).name fmap(file).magnitude];
		bids(subj).fmap(file).phasediff=[bids_dir bids(subj).name '/fmap/' bids(subj).name fmap(file).phasediff];
	end

	for file=1:length(func)
		bids(subj).func(file).fname=[bids_dir bids(subj).name '/func/' bids(subj).name func(file).fname];
	end

end

% load in information from .json
for subj=1:num_subj

	for file=1:length(fmap)
		if exist([bids_dir bids(subj).name '/fmap/' bids(subj).name fmap(file).phasediff '.json'],'file')
			params=json2mat([bids_dir bids(subj).name '/fmap/' bids(subj).name fmap(file).phasediff '.json']);
		else
			params=json2mat([bids_dir fmap(file).phasediff(2:end) '.json']);
		end
		bids(subj).fmap(file).params=params;
	end

	for file=1:length(func)
		if exist([bids_dir bids(subj).name '/func/' bids(subj).name func(file).fname '.json'],'file')
			params=json2mat([bids_dir bids(subj).name '/func/' bids(subj).name func(file).fname '.json']);
		else
			params=json2mat([bids_dir func(file).fname(2:end) '.json']);
		end
		bids(subj).func(file).params=params;
	end
	
end








