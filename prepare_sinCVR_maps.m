function bids=prepare_sinCVR_maps(bids,subj)
	
	%using full data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(1).feat 'stats/sigmasquareds']);
		
	%read in resp data and estimate change in PetCO2
	events_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_events.tsv'];
	fid=fopen(events_in,'r');
	events=textscan(fid,'%f %f %s','delimiter','\t','headerlines',1);
	fclose(fid);
	
	row=find(strcmp(events{3},'sinCVR'));
	bbb_in=[bids(subj).dir bids(subj).name '/resp/' bids(subj).name '_respdata_bbb.tsv'];
	bbb=importdata(bbb_in,'\t');
	ind1=find(bbb.data(:,1)>events{1}(row),1,'first');
	ind2=find(bbb.data(:,1)<events{2}(row),1,'last');
	time=bbb.data(ind1:ind2,1)-bbb.data(ind1,1);
	time=time.*60;
	sinpetco2=bbb.data(ind1:ind2,2);	
	timei=(1:210)'.*2;	
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;
	
	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');
	
	bids(subj).func(1).results(1).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(1).petco2basevar=sinstd(3);
	bids(subj).func(1).results(1).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(1).petco2deltavar=sinvardeltapetco2;
			
	%estimate magnitude and errors
	sincvr_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr_mag(isnan(sincvr_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr_magvar(isnan(sincvr_magvar))=0;
	
	sincvr_magrsd=sqrt(sincvr_magvar)./sincvr_mag;
	
	%estimate phase and errors
	sincvr_pha=atan2(cope1,cope2); %use atan2 for more robustness
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr_phavar=varcopediv./(1+(cope1./cope2).^2).^2;
	sincvr_phavar(isnan(sincvr_phavar))=0;
		
	%write out images
	bids(subj).func(1).results(1).name='cvr magnitude/phase maps from sinCVR - full data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(1).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-mag'];
	save_avw(sincvr_mag,bids(subj).func(1).results(1).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(1).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-magvar'];
	save_avw(sincvr_magvar,bids(subj).func(1).results(1).cvr_magvar,'f',scales');

	bids(subj).func(1).results(1).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-magrsd'];
	save_avw(sincvr_magrsd,bids(subj).func(1).results(1).cvr_magrsd,'f',scales');
	
	bids(subj).func(1).results(1).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-pha'];
	save_avw(sincvr_pha,bids(subj).func(1).results(1).cvr_pha,'f',scales');

	bids(subj).func(1).results(1).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-phavar'];
	save_avw(sincvr_phavar,bids(subj).func(1).results(1).cvr_phavar,'f',scales');
		
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(1).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(1).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(1).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(1).global_pha=global_pha;
	
	sincvr_phan=sincvr_pha;
	sincvr_phan(sincvr_phan~=0)=sincvr_phan(sincvr_phan~=0)-global_pha;
	sincvr_phan(sincvr_phan<-pi)=sincvr_phan(sincvr_phan<-pi)+2*pi;
	sincvr_phan(sincvr_phan>pi)=sincvr_phan(sincvr_phan>pi)-2*pi;
	
	save_avw(sincvr_phan,[bids(subj).func(1).results(1).cvr_pha 'n'],'f',scales');
	
	
	%using 6 mins data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(2).feat 'stats/sigmasquareds']);

	%estimate change in PetCO2	
	timei=(1:180)'.*2;
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;
	
	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');
	
	bids(subj).func(1).results(2).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(2).petco2basevar=sinstd(3);
	bids(subj).func(1).results(2).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(2).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr6min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr6min_mag(isnan(sincvr6min_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr6min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr6min_magvar(isnan(sincvr6min_magvar))=0;
	
	sincvr6min_magrsd=sqrt(sincvr6min_magvar)./sincvr6min_mag;

	%estimate phase and errors
	sincvr6min_pha=atan2(cope1,cope2);	
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr6min_phavar=varcopediv./(1+(cope1./cope2).^2).^2;	
	sincvr6min_phavar(isnan(sincvr6min_phavar))=0;
	
	%write out images
	bids(subj).func(1).results(2).name='cvr magnitude/phase maps from sinCVR - 6 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(2).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-mag'];
	save_avw(sincvr6min_mag,bids(subj).func(1).results(2).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(2).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-magvar'];
	save_avw(sincvr6min_magvar,bids(subj).func(1).results(2).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(2).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-magrsd'];
	save_avw(sincvr6min_magrsd,bids(subj).func(1).results(2).cvr_magrsd,'f',scales');
		
	bids(subj).func(1).results(2).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-pha'];
	save_avw(sincvr6min_pha,bids(subj).func(1).results(2).cvr_pha,'f',scales');
	
	bids(subj).func(1).results(2).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-phavar'];
	save_avw(sincvr6min_phavar,bids(subj).func(1).results(2).cvr_phavar,'f',scales');
		
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(2).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(2).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(2).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(2).global_pha=global_pha;
	
	sincvr6min_phan=sincvr6min_pha;
	sincvr6min_phan(sincvr6min_phan~=0)=sincvr6min_phan(sincvr6min_phan~=0)-global_pha;
	sincvr6min_phan(sincvr6min_phan<-pi)=sincvr6min_phan(sincvr6min_phan<-pi)+2*pi;
	sincvr6min_phan(sincvr6min_phan>pi)=sincvr6min_phan(sincvr6min_phan>pi)+2*pi;
	
	save_avw(sincvr6min_phan,[bids(subj).func(1).results(2).cvr_pha 'n'],'f',scales');
	
	
	%using 5 mins data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(3).feat 'stats/sigmasquareds']);

	%estimate change in PetCO2	
	timei=(1:150)'.*2;
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;
	
	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');
	
	bids(subj).func(1).results(3).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(3).petco2basevar=sinstd(3);
	bids(subj).func(1).results(3).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(3).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr5min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr5min_mag(isnan(sincvr5min_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr5min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr5min_magvar(isnan(sincvr5min_magvar))=0;
	
	sincvr5min_magrsd=sqrt(sincvr5min_magvar)./sincvr5min_mag;
	
	%estimate phase and errors
	sincvr5min_pha=atan2(cope1,cope2);	
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr5min_phavar=varcopediv./(1+(cope1./cope2).^2).^2;	
	sincvr5min_phavar(isnan(sincvr5min_phavar))=0;
	
	%write out images
	bids(subj).func(1).results(3).name='cvr magnitude/phase maps from sinCVR - 5 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(3).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-mag'];
	save_avw(sincvr5min_mag,bids(subj).func(1).results(3).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(3).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-magvar'];
	save_avw(sincvr5min_magvar,bids(subj).func(1).results(3).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(3).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-magrsd'];
	save_avw(sincvr5min_magrsd,bids(subj).func(1).results(3).cvr_magrsd,'f',scales');
		
	bids(subj).func(1).results(3).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-pha'];
	save_avw(sincvr5min_pha,bids(subj).func(1).results(3).cvr_pha,'f',scales');
	
	bids(subj).func(1).results(3).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-phavar'];
	save_avw(sincvr5min_phavar,bids(subj).func(1).results(3).cvr_phavar,'f',scales');
		
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(3).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(3).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(3).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(3).global_pha=global_pha;
	
	sincvr5min_phan=sincvr5min_pha;
	sincvr5min_phan(sincvr5min_phan~=0)=sincvr5min_phan(sincvr5min_phan~=0)-global_pha;
	sincvr5min_phan(sincvr5min_phan<-pi)=sincvr5min_phan(sincvr5min_phan<-pi)+2*pi;
	sincvr5min_phan(sincvr5min_phan>pi)=sincvr5min_phan(sincvr5min_phan>pi)+2*pi;
	
	save_avw(sincvr5min_phan,[bids(subj).func(1).results(3).cvr_pha 'n'],'f',scales');
	
	
	%using 4 mins data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(4).feat 'stats/sigmasquareds']);

	%estimate change in PetCO2	
	timei=(1:120)'.*2;
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;
	
	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');
	
	bids(subj).func(1).results(4).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(4).petco2basevar=sinstd(3);
	bids(subj).func(1).results(4).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(4).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr4min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr4min_mag(isnan(sincvr4min_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr4min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr4min_magvar(isnan(sincvr4min_magvar))=0;
	
	sincvr4min_magrsd=sqrt(sincvr4min_magvar)./sincvr4min_mag;
	
	%estimate phase and errors
	sincvr4min_pha=atan2(cope1,cope2);	
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr4min_phavar=varcopediv./(1+(cope1./cope2).^2).^2;	
	sincvr4min_phavar(isnan(sincvr4min_phavar))=0;
	
	%write out images
	bids(subj).func(1).results(4).name='cvr magnitude/phase maps from sinCVR - 5 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(4).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-mag'];
	save_avw(sincvr4min_mag,bids(subj).func(1).results(4).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(4).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-magvar'];
	save_avw(sincvr4min_magvar,bids(subj).func(1).results(4).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(4).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-magrsd'];
	save_avw(sincvr4min_magrsd,bids(subj).func(1).results(4).cvr_magrsd,'f',scales');
		
	bids(subj).func(1).results(4).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-pha'];
	save_avw(sincvr4min_pha,bids(subj).func(1).results(4).cvr_pha,'f',scales');
	
	bids(subj).func(1).results(4).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-phavar'];
	save_avw(sincvr4min_phavar,bids(subj).func(1).results(4).cvr_phavar,'f',scales');
		
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(4).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(4).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(4).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(4).global_pha=global_pha;
	
	sincvr4min_phan=sincvr4min_pha;
	sincvr4min_phan(sincvr4min_phan~=0)=sincvr4min_phan(sincvr4min_phan~=0)-global_pha;
	sincvr4min_phan(sincvr4min_phan<-pi)=sincvr4min_phan(sincvr4min_phan<-pi)+2*pi;
	sincvr4min_phan(sincvr4min_phan>pi)=sincvr4min_phan(sincvr4min_phan>pi)+2*pi;
	
	save_avw(sincvr4min_phan,[bids(subj).func(1).results(4).cvr_pha 'n'],'f',scales');
	
		
	%using 3 mins data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(5).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(5).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(5).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(5).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(5).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(5).feat 'stats/sigmasquareds']);
	
	%estimate change in PetCO2
	timei=(1:90)'.*2;
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;

	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');

	bids(subj).func(1).results(5).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(5).petco2basevar=sinstd(3);
	bids(subj).func(1).results(5).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(5).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr3min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr3min_mag(isnan(sincvr3min_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr3min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr3min_magvar(isnan(sincvr3min_magvar))=0;
	
	sincvr3min_magrsd=sqrt(sincvr3min_magvar)./sincvr3min_mag;
	
	%estimate phase and errors
	sincvr3min_pha=atan2(cope1,cope2);
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr3min_phavar=varcopediv./(1+(cope1./cope2).^2).^2;	
	sincvr3min_phavar(isnan(sincvr3min_phavar))=0;
	
	%write out images 
	bids(subj).func(1).results(5).name='cvr magnitude/phase maps from sinCVR - 3 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(5).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-mag'];
	save_avw(sincvr3min_mag,bids(subj).func(1).results(5).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(5).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-magvar'];
	save_avw(sincvr3min_magvar,bids(subj).func(1).results(5).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(5).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-magrsd'];
	save_avw(sincvr3min_magrsd,bids(subj).func(1).results(5).cvr_magrsd,'f',scales');
		
	bids(subj).func(1).results(5).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-pha'];
	save_avw(sincvr3min_pha,bids(subj).func(1).results(5).cvr_pha,'f',scales');
	
	bids(subj).func(1).results(5).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-phavar'];
	save_avw(sincvr3min_phavar,bids(subj).func(1).results(5).cvr_phavar,'f',scales');	
	
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(5).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(5).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(5).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(5).global_pha=global_pha;
	
	sincvr3min_phan=sincvr3min_pha;
	sincvr3min_phan(sincvr3min_phan~=0)=sincvr3min_phan(sincvr3min_phan~=0)-global_pha;
	sincvr3min_phan(sincvr3min_phan<-pi)=sincvr3min_phan(sincvr3min_phan<-pi)+2*pi;
	sincvr3min_phan(sincvr3min_phan>pi)=sincvr3min_phan(sincvr3min_phan>pi)+2*pi;
	
	save_avw(sincvr3min_phan,[bids(subj).func(1).results(5).cvr_pha 'n'],'f',scales');
	
		
	%using 2 mins data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(6).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(6).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(6).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(6).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(6).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(6).feat 'stats/sigmasquareds']);
	
	%estimate change in PetCO2
	timei=(1:60)'.*2;
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;

	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');

	bids(subj).func(1).results(6).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(6).petco2basevar=sinstd(3);
	bids(subj).func(1).results(6).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(6).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr2min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr2min_mag(isnan(sincvr2min_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr2min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr2min_magvar(isnan(sincvr2min_magvar))=0;
	
	sincvr2min_magrsd=sqrt(sincvr2min_magvar)./sincvr2min_mag;
	
	%estimate phase and errors
	sincvr2min_pha=atan2(cope1,cope2);
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr2min_phavar=varcopediv./(1+(cope1./cope2).^2).^2;	
	sincvr2min_phavar(isnan(sincvr2min_phavar))=0;
	
	%write out images 
	bids(subj).func(1).results(6).name='cvr magnitude/phase maps from sinCVR - 2 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(6).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-mag'];
	save_avw(sincvr2min_mag,bids(subj).func(1).results(6).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(6).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-magvar'];
	save_avw(sincvr2min_magvar,bids(subj).func(1).results(6).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(6).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-magrsd'];
	save_avw(sincvr2min_magrsd,bids(subj).func(1).results(6).cvr_magrsd,'f',scales');
		
	bids(subj).func(1).results(6).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-pha'];
	save_avw(sincvr2min_pha,bids(subj).func(1).results(6).cvr_pha,'f',scales');
	
	bids(subj).func(1).results(6).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-phavar'];
	save_avw(sincvr2min_phavar,bids(subj).func(1).results(6).cvr_phavar,'f',scales');	
	
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(6).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(6).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(6).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(6).global_pha=global_pha;
	
	sincvr2min_phan=sincvr2min_pha;
	sincvr2min_phan(sincvr2min_phan~=0)=sincvr2min_phan(sincvr2min_phan~=0)-global_pha;
	sincvr2min_phan(sincvr2min_phan<-pi)=sincvr2min_phan(sincvr2min_phan<-pi)+2*pi;
	sincvr2min_phan(sincvr2min_phan>pi)=sincvr2min_phan(sincvr2min_phan>pi)+2*pi;
	
	save_avw(sincvr2min_phan,[bids(subj).func(1).results(6).cvr_pha 'n'],'f',scales');
		
		
	%using 1 min data
	
	%read in source images
	[cope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(7).feat 'stats/cope1']);
	[cope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(7).feat 'stats/cope2']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(7).feat 'mean_func']);
	[varcope1 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(7).feat 'stats/varcope1']);
	[varcope2 dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(7).feat 'stats/varcope2']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(7).feat 'stats/sigmasquareds']);
	
	%estimate change in PetCO2
	timei=(1:30)'.*2;
	sinpetco2i=interp1(time,sinpetco2,timei,'linear','extrap');
	sinstim=[sin(timei./60.*2.*pi) cos(timei./60.*2.*pi)];
	sinX=[sinstim ones(size(timei))];
	[sina sinstd]=lscov(sinX,sinpetco2i);
	sindeltapetco2=sqrt(sum(sina(1:2).^2)).*2;
	sinvardeltapetco2=(sina(1)^2./(sina(1)^2+sina(2)^2))*sinstd(1).^2+(sina(2)^2./(sina(1)^2+sina(2)^2))*sinstd(2).^2;

	%figure;
	%plot(sinpetco2i,'x');
	%hold on;
	%plot(sinX*sina,'-');

	bids(subj).func(1).results(7).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(7).petco2basevar=sinstd(3);
	bids(subj).func(1).results(7).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(7).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr1min_mag=sqrt(cope1.^2+cope2.^2)./mean_func;
	sincvr1min_mag(isnan(sincvr1min_mag))=0;
	
	cope=sqrt(cope1.^2+cope2.^2);
	varcope=(cope1.^2./(cope1.^2+cope2.^2)).*varcope1+(cope2.^2./(cope1.^2+cope2.^2)).*varcope2;
	sincvr1min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr1min_magvar(isnan(sincvr1min_magvar))=0;
	
	sincvr1min_magrsd=sqrt(sincvr1min_magvar)./sincvr1min_mag;
	
	%estimate phase and errors
	sincvr1min_pha=atan2(cope1,cope2);
	
	varcopediv=(cope1./cope2).^2.*(varcope1./cope1.^2+varcope2./cope2.^2);
	sincvr1min_phavar=varcopediv./(1+(cope1./cope2).^2).^2;	
	sincvr1min_phavar(isnan(sincvr1min_phavar))=0;
	
	%write out images 
	bids(subj).func(1).results(7).name='cvr magnitude/phase maps from sinCVR - 1 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(7).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-mag'];
	save_avw(sincvr1min_mag,bids(subj).func(1).results(7).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(7).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-magvar'];
	save_avw(sincvr1min_magvar,bids(subj).func(1).results(7).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(7).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-magrsd'];
	save_avw(sincvr1min_magrsd,bids(subj).func(1).results(7).cvr_magrsd,'f',scales');
		
	bids(subj).func(1).results(7).cvr_pha=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-pha'];
	save_avw(sincvr1min_pha,bids(subj).func(1).results(7).cvr_pha,'f',scales');
	
	bids(subj).func(1).results(7).cvr_phavar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-phavar'];
	save_avw(sincvr1min_phavar,bids(subj).func(1).results(7).cvr_phavar,'f',scales');	
	
	%adjust phase values to global phase across brain and unwrap to be comparable with Toronto protocol
	[func fdims]=read_avw([bids(subj).func(1).analysis(7).feat 'filtered_func_data.nii.gz']);
	funcr=reshape(func,prod(fdims(1:3)),fdims(4));
	mask=read_avw([bids(subj).func(1).analysis(7).feat 'mask']);
	maskr=mask(:);                                                               
	tc=mean(funcr(mask>0,:));
	tcm=tc-mean(tc);
	ev=dlmread([bids(subj).func(1).analysis(7).feat 'design.mat'],'\t',5,0);
	X=[ev(:,1:2)];
	a=X\tcm';
	global_pha=atan2(a(1),a(2));
	bids(subj).func(1).results(7).global_pha=global_pha;
	
	sincvr1min_phan=sincvr1min_pha;
	sincvr1min_phan(sincvr1min_phan~=0)=sincvr1min_phan(sincvr1min_phan~=0)-global_pha;
	sincvr1min_phan(sincvr1min_phan<-pi)=sincvr1min_phan(sincvr1min_phan<-pi)+2*pi;
	sincvr1min_phan(sincvr1min_phan>pi)=sincvr1min_phan(sincvr1min_phan>pi)+2*pi;
	
	save_avw(sincvr1min_phan,[bids(subj).func(1).results(7).cvr_pha 'n'],'f',scales');