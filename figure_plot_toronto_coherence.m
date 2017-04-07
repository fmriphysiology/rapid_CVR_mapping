function figure_plot_toronto_coherence(bids)
	
	for subj=1:9
	
		[C(:,subj) F]=plot_toronto_coherence(bids,subj);
		
	end
	
	figure;
	plot(F,mean(C,2),'k');
	hold on;
	plot(F,mean(C,2)+std(C,[],2),'k--');
	plot(F,mean(C,2)-std(C,[],2),'k--');
	axis square;
	title('Whole brain coherence as a function of frequency');
	
	
	%keyboard;