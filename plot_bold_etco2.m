function [gmout wmout timei petco2out timei2]=plot_bold_etco2(bids,subj);
	
	%register gm/wm masks
	in=bids(subj).anat.gm;
	ref=[bids(subj).func(1).analysis(1).feat 'mean_func'];
	out=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-sinusoidCVR_bold_gm'];
	xfm=[bids(subj).func(1).analysis(1).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -init ' xfm ' -applyxfm']);

	in=bids(subj).anat.wm;
	ref=[bids(subj).func(1).analysis(1).feat 'mean_func'];
	out=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-sinusoidCVR_bold_wm'];
	xfm=[bids(subj).func(1).analysis(1).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -init ' xfm ' -applyxfm']);	

	in=bids(subj).anat.gm;
	ref=[bids(subj).func(2).analysis(2).feat 'mean_func'];
	out=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-torontoCVR_bold_gm'];
	xfm=[bids(subj).func(2).analysis(2).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -init ' xfm ' -applyxfm']);

	in=bids(subj).anat.wm;
	ref=[bids(subj).func(2).analysis(2).feat 'mean_func'];
	out=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-torontoCVR_bold_wm'];
	xfm=[bids(subj).func(2).analysis(2).feat 'reg/highres2example_func.mat'];
	
	status=system(['flirt -in ' in ' -ref ' ref ' -out ' out ' -init ' xfm ' -applyxfm']);	
	
	%load in sinusoid images
	img=read_avw(bids(subj).func(1).fname);
	gm=read_avw([bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-sinusoidCVR_bold_gm']);
	wm=read_avw([bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-sinusoidCVR_bold_wm']);	
	
	events_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_events.tsv'];
	fid=fopen(events_in,'r');
	events=textscan(fid,'%f %f %s','delimiter','\t','headerlines',1);
	fclose(fid);
	row=find(strcmp(events{3},'sinCVR'));
	bbb_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_bbb.tsv'];
	bbb=importdata(bbb_in,'\t');
	ind1=find(bbb.data(:,1)>events{1}(row),1,'first');
	ind2=find(bbb.data(:,1)<events{2}(row),1,'last');
	ind1a=find(bbb.data(:,1)>(events{1}(row)-1),1,'first');
	ind2a=find(bbb.data(:,1)<(events{2}(row)+1),1,'last');
	time=bbb.data(ind1a:ind2a,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1a:ind2a,2);
	timei=(1:210)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	timei2=(-30:240)'.*2;
	petco2i2=interp1(time,petco2,timei2,'linear','extrap');
	
	ind1a=find(bbb.data(:,1)>(events{1}(row)-1),1,'first');
	ind2a=find(bbb.data(:,1)<(events{2}(row)+1),1,'last');
	
	imgr=reshape(img,64*64*24,210);
	gmr=gm(:);
	wmr=wm(:);
	
	gmtc=mean(imgr(gmr>0.5,:))./mean(mean(imgr(gmr>0.5,:))).*100-100;
	wmtc=mean(imgr(wmr>0.5,:))./mean(mean(imgr(wmr>0.5,:))).*100-100;

	gmout(:,1)=gmtc;
	wmout(:,1)=wmtc;
	petco2out(:,1)=petco2i2;

	%load in toronto images
	img=read_avw(bids(subj).func(2).fname);
	gm=read_avw([bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-torontoCVR_bold_gm']);
	wm=read_avw([bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).name '_task-torontoCVR_bold_wm']);	
	
	events_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_events.tsv'];
	fid=fopen(events_in,'r');
	events=textscan(fid,'%f %f %s','delimiter','\t','headerlines',1);
	fclose(fid);
	row=find(strcmp(events{3},'torontoCVR'));
	bbb_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_bbb.tsv'];
	bbb=importdata(bbb_in,'\t');
	ind1=find(bbb.data(:,1)>events{1}(row),1,'first');
	ind2=find(bbb.data(:,1)<events{2}(row),1,'last');
	ind1a=find(bbb.data(:,1)>(events{1}(row)-1),1,'first');
	ind2a=find(bbb.data(:,1)<(events{2}(row)+1),1,'last');
	time=bbb.data(ind1a:ind2a,1)-bbb.data(ind1,1);
	time=time.*60;
	petco2=bbb.data(ind1a:ind2a,2);
	timei=(1:210)'.*2;
	petco2i=interp1(time,petco2,timei,'linear','extrap');
	timei2=(-30:240)'.*2;
	petco2i2=interp1(time,petco2,timei2,'linear','extrap');
	
	ind1a=find(bbb.data(:,1)>(events{1}(row)-1),1,'first');
	ind2a=find(bbb.data(:,1)<(events{2}(row)+1),1,'last');
	
	imgr=reshape(img,64*64*24,210);
	gmr=gm(:);
	wmr=wm(:);
	
	gmtc=mean(imgr(gmr>0.5,:))./mean(mean(imgr(gmr>0.5,:))).*100-100;
	wmtc=mean(imgr(wmr>0.5,:))./mean(mean(imgr(wmr>0.5,:))).*100-100;

	gmout(:,2)=gmtc;
	wmout(:,2)=wmtc;
	petco2out(:,2)=petco2i2;