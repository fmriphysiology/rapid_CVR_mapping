function [C F phi]=plot_toronto_coherence(bids,subj);
	
	%cvr phase from transfer function analysis (no errors since no defined method for calculating variance)
	ev=dlmread([bids(subj).func(2).analysis(2).feat 'design.mat'],'\t',5,0);
	[func_data func_dims]=read_avw([bids(subj).func(2).analysis(2).feat 'filtered_func_data']);
	func=reshape(func_data,prod(func_dims(1:3)),func_dims(4));
	[mask_data]=read_avw([bids(subj).func(2).analysis(2).feat 'mask']);
	mask=mask_data(:);
	braintc=mean(func(mask>0,:));

	%[H F]=tfestimate(ev(:,1)-mean(ev(:,1)),func-repmat(mean(func,1),func_dims(4),1),50,[],[],1/2);
	
	[pevbold f]=cpsd(ev(:,1)-mean(ev(:,1)),braintc-mean(braintc),50,[],[],1/2);
	[pev f]=pwelch(ev(:,1)-mean(ev(:,1)),50,[],[],1/2);
	H=pevbold./pev;
	
	%freq=find(f>0.01,1,'first');
	%freq=find(f>1/60,1,'first');
	%torcvr_pha=reshape(atan2(imag(H(freq,:)),real(H(freq,:))),func_dims(1),func_dims(2),func_dims(3)); %use atan2 for more robustness

	freq=find(f>0.01,1,'first');
	phi(:,1)=atan2(imag(H(freq,:)),real(H(freq,:)));
	freq=find(f>1/60,1,'first');
	phi(:,2)=atan2(imag(H(freq,:)),real(H(freq,:)));

	[C F]=mscohere(ev(:,1)-mean(ev(:,1)),braintc-mean(braintc),50,[],[],1/2);
	
	%keyboard;
	