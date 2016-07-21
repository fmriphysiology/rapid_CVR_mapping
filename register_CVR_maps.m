function bids=register_CVR_maps(bids,subj)

	%register all mean_func to toronto mean_func

	%register sinusoid to toronto
	in=[bids(subj).func(1).analysis(1).feat 'mean_func'];
	ref=[bids(subj).func(2).analysis(2).feat 'mean_func'];
	omat{1}=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' 'sin2tor.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -dof 6 -omat ' omat{1}]);
	
	%register sinusoid 5min to toronto
	in=[bids(subj).func(1).analysis(2).feat 'mean_func'];
	ref=[bids(subj).func(2).analysis(2).feat 'mean_func'];
	omat{2}=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' 'sin5min2tor.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -dof 6 -omat ' omat{2}]);
		
	%register sinusoid 3min to toronto
	in=[bids(subj).func(1).analysis(3).feat 'mean_func'];
	ref=[bids(subj).func(2).analysis(2).feat 'mean_func'];
	omat{3}=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' 'sin3min2tor.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -dof 6 -omat ' omat{3}]);
	
	%transform all cvr mag maps to toronto space
	for k=1:3
		in=bids(subj).func(1).results(k).cvr_mag;
		out=[bids(subj).func(1).results(k).cvr_mag '_reg'];
		ref=bids(subj).func(2).results(1).cvr_mag;
	
		status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat{k}]);
	end
	
	%transform all cvr mag variance maps to toronto space
	for k=1:3
		in=bids(subj).func(1).results(k).cvr_magvar;
		out=[bids(subj).func(1).results(k).cvr_magvar '_reg'];
		ref=bids(subj).func(2).results(1).cvr_mag;
	
		status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat{k}]);
	end

	%transform all cvr phase maps to toronto
	for k=1:3
		in=bids(subj).func(1).results(k).cvr_pha;
		out=[bids(subj).func(1).results(k).cvr_pha '_reg'];
		ref=bids(subj).func(2).results(1).cvr_pha;
	
		status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat{k}]);
	end
	
	%transform all cvr phase (adjusted) maps to toronto
	for k=1:3
		in=[bids(subj).func(1).results(k).cvr_pha 'n'];
		out=[bids(subj).func(1).results(k).cvr_pha 'n_reg'];
		ref=bids(subj).func(2).results(1).cvr_pha;
	
		status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat{k}]);
	end
	
	%transform all cvr phase variance maps to toronto
	for k=1:3
		in=[bids(subj).func(1).results(k).cvr_phavar];
		out=[bids(subj).func(1).results(k).cvr_phavar '_reg'];
		ref=bids(subj).func(2).results(1).cvr_pha;
	
		status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat{k}]);
	end	
	
	%transform gm, wm and csf pve masks to toronto space
	in=bids(subj).anat.gm;
	out=[bids(subj).anat.gm '_reg'];
	ref=bids(subj).func(2).results(1).cvr_mag;
	omat=[bids(subj).func(1).analysis(2).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat]);
	
	in=bids(subj).anat.wm;
	s=regexp(in,'/');
	out=[bids(subj).anat.wm '_reg'];
	ref=bids(subj).func(2).results(1).cvr_mag;
	omat=[bids(subj).func(1).analysis(2).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat]);
	
	in=bids(subj).anat.csf;
	s=regexp(in,'/');
	out=[bids(subj).anat.csf '_reg'];
	ref=bids(subj).func(2).results(1).cvr_mag;
	omat=[bids(subj).func(1).analysis(2).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' omat]);
	
		
	