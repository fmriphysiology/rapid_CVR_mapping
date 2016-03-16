function bids=prepare_anat(bids,subj);

	status=system([bids(subj).dir 'derivatives/code/prepare_anat.sh ' bids(subj).dir ' ' bids(subj).name]);
	if status==0
		bids(subj).anat.bet=[bids(subj).dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain'];
		bids(subj).anat.csf=[bids(subj).dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain_pve_0'];
		bids(subj).anat.gm=[bids(subj).dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain_pve_1'];
		bids(subj).anat.wm=[bids(subj).dir 'derivatives/' bids(subj).name '/anat/' bids(subj).name '_T1w_brain_pve_2'];
	end