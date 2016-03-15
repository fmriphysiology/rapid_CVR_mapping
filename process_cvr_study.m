parfor subj=1:10;

	bids_dir='/Users/nickb/Analysis/fmrib/cvr_study/';

	% initialise bids data structure
	bids=bids_init(bids_dir);

	% prepare fieldmaps
	bids=prepare_fieldmap(bids,subj);

	% prepare anatomical images
	bids=prepare_anat(bids,subj);

	% prepare slice order file
	bids=prepare_sliceorder(bids,subj);

	% prepare and run 1st pass sinCVR feat design
	bids=process_sinCVR_firstpass(bids,subj);

	% prepare 1st pass torontoCVR feat design
	bids=process_torontoCVR_firstpass(bids,subj);
	
	% prepare 2nd pass sinCVR feat designs - data reduction
	bids=process_sinCVR_secondpass(bids,subj);

	% prepare 2nd pass torontoCVR feat designs - subject specific PetCO2 regressor
	bids=process_torontoCVR_secondpass(bids,subj);
	
	% prepare 3rd pass torontoCVR transfer function analysis???
	
end