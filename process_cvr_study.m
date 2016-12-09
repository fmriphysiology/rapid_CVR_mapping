%add dependent paths
addpath ~/Documents/MATLAB/json4mat
addpath /usr/local/fsl/etc/matlab

bids_dir='/Users/nickb/Analysis/fmrib/cvr_study/';

% initialise bids data structure
bids=bids_init(bids_dir);

for subj=1:10;

	disp(['Processing subject ' num2str(subj) '...']);

	% create output file structure
	bids=create_output_dirs(bids,subj);

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
	
	% prepare sinCVR maps
	bids=prepare_sinCVR_maps(bids,subj);
	
	% prepare torontoCVR_maps(bids,subj);
	bids=prepare_torontoCVR_maps(bids,subj);
	
	%register CVR maps and tissue masks to Toronto functional space
	bids=register_CVR_maps(bids,subj);
	
	% save progress
	bids=save_progress(bids,subj);
	
end