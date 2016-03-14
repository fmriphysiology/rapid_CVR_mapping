function bids=process_fieldmap(bids,subj)

	deltaTE=(str2num(bids(subj).fmap.params.EchoTime2)-str2num(bids(subj).fmap.params.EchoTime1))*1000;
	status=system([bids(subj).dir 'derivatives/code/prep_fmap.sh ' bids(subj).dir ' ' bids(subj).name ' ' num2str(deltaTE)]);
	if status==0
		bids(subj).fmap.fieldmap=[bids(subj).dir 'derivatives/' bids(subj).name '/fmap/' bids(subj).name '_fieldmap'];
		bids(subj).fmap.fieldmap_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/fmap/' bids(subj).name '_magnitude_e2_brain_ero'];
	end