subj=5;

bids_dir='/Users/nickb/Analysis/fmrib/cvr_study/';

% initialise bids data structure
bids=bids_init(bids_dir);

% prepare fieldmaps
deltaTE=(str2num(bids(subj).fmap(file).params.EchoTime2)-str2num(bids(subj).fmap.params.EchoTime1))*1000;
status=system([bids_dir 'derivatives/code/prep_fmap.sh ' bids_dir ' ' bids(subj).name ' ' num2str(deltaTE)]);
if status==0
	bids(subj).fmap.fieldmap=[bids_dir 'derivatives/' bids(subj).name '/fmap/' bids(subj).name '_fieldmap'];
end

% prepare anatomical images
status=system([bids_dir 'derivatives/code/prep_anat.sh ' bids_dir ' ' bids(subj).name]);
if status==0
	bids(subj).anat.bet=[bids_dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain'];
end

% prepare slice timing file
for file=1:length(bids(subj).func)
	[timing sliceorder]=sort(bids(subj).func(file).params.SliceTiming');
	if ~exist([bids_dir 'derivatives/' bids(subj).name])
		system(['mkdir ' bids_dir 'derivatives/' bids(subj).name]);
	end
	if ~exist([bids_dir 'derivatives/' bids(subj).name '/func'])
		system(['mkdir ' bids_dir 'derivatives/' bids(subj).name '/func']);
	end
	s=regexp(bids(subj).func(file).fname,'/');
	dlmwrite([bids_dir 'derivatives/' bids(subj).name '/func/' bids(1).func(1).fname(s(end)+1:end) '_sliceorder.tsv'],sliceorder);
	if exist([bids_dir 'derivatives/' bids(subj).name '/func/' bids(1).func(1).fname(s(end)+1:end) '_sliceorder.tsv'],'file');
		bids(subj).func(file).sliceorder=[bids_dir 'derivatives/' bids(subj).name '/func/' bids(1).func(1).fname(s(end)+1:end) '_sliceorder.tsv'];
	end
end


