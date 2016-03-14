function bids=prepare_anat(bids,subj);

	status=system([bids(subj).dir 'derivatives/code/prep_anat.sh ' bids(subj).dir ' ' bids(subj).name]);
	if status==0
		bids(subj).anat.bet=[bids(subj).dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain'];
	end