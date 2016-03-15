function bids=save_progress(bids,subj);

	save([bids(subj).dir 'derivatives/bids.mat'],'bids');
	