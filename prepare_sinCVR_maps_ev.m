function bids=prepare_sinCVR_maps_ev(bids,subj)
	
	%using full data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(8).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(8).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(8).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(8).feat 'stats/sigmasquareds']);
		
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
	
	bids(subj).func(1).results(8).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(8).petco2basevar=sinstd(3);
	bids(subj).func(1).results(8).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(8).petco2deltavar=sinvardeltapetco2;
			
	%estimate magnitude and errors
	sincvr_mag=cope./mean_func;
	sincvr_mag(isnan(sincvr_mag))=0;
	
	sincvr_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr_magvar(isnan(sincvr_magvar))=0;
	
	sincvr_magrsd=sqrt(sincvr_magvar)./sincvr_mag;
		
	%write out images
	bids(subj).func(1).results(8).name='cvr magnitude/phase maps from sinCVR - full data ev';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(8).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-ev-mag'];
	save_avw(sincvr_mag,bids(subj).func(1).results(8).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(8).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-ev-magvar'];
	save_avw(sincvr_magvar,bids(subj).func(1).results(8).cvr_magvar,'f',scales');

	bids(subj).func(1).results(8).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr-ev-magrsd'];
	save_avw(sincvr_magrsd,bids(subj).func(1).results(8).cvr_magrsd,'f',scales');
	
	%using 6 mins data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(9).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(9).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(9).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(9).feat 'stats/sigmasquareds']);

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
	
	bids(subj).func(1).results(9).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(9).petco2basevar=sinstd(3);
	bids(subj).func(1).results(9).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(9).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr6min_mag=cope./mean_func;
	sincvr6min_mag(isnan(sincvr6min_mag))=0;
	
	sincvr6min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr6min_magvar(isnan(sincvr6min_magvar))=0;
	
	sincvr6min_magrsd=sqrt(sincvr6min_magvar)./sincvr6min_mag;
	
	%write out images
	bids(subj).func(1).results(9).name='cvr magnitude/phase maps from sinCVR - 6 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(9).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-ev-mag'];
	save_avw(sincvr6min_mag,bids(subj).func(1).results(9).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(9).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-ev-magvar'];
	save_avw(sincvr6min_magvar,bids(subj).func(1).results(9).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(9).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr6min-ev-magrsd'];
	save_avw(sincvr6min_magrsd,bids(subj).func(1).results(9).cvr_magrsd,'f',scales');
	
	%using 5 mins data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(10).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(10).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(10).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(10).feat 'stats/sigmasquareds']);

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
	
	bids(subj).func(1).results(10).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(10).petco2basevar=sinstd(3);
	bids(subj).func(1).results(10).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(10).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr5min_mag=cope./mean_func;
	sincvr5min_mag(isnan(sincvr5min_mag))=0;
	
	sincvr5min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr5min_magvar(isnan(sincvr5min_magvar))=0;
	
	sincvr5min_magrsd=sqrt(sincvr5min_magvar)./sincvr5min_mag;
	
	%write out images
	bids(subj).func(1).results(10).name='cvr magnitude/phase maps from sinCVR - 5 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(10).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-ev-mag'];
	save_avw(sincvr5min_mag,bids(subj).func(1).results(10).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(10).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-ev-magvar'];
	save_avw(sincvr5min_magvar,bids(subj).func(1).results(10).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(10).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr5min-ev-magrsd'];
	save_avw(sincvr5min_magrsd,bids(subj).func(1).results(10).cvr_magrsd,'f',scales');

	%using 4 mins data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(11).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(11).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(11).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(11).feat 'stats/sigmasquareds']);

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
	
	bids(subj).func(1).results(11).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(11).petco2basevar=sinstd(3);
	bids(subj).func(1).results(11).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(11).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr4min_mag=cope./mean_func;
	sincvr4min_mag(isnan(sincvr4min_mag))=0;
	
	sincvr4min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr4min_magvar(isnan(sincvr4min_magvar))=0;
	
	sincvr4min_magrsd=sqrt(sincvr4min_magvar)./sincvr4min_mag;
	
	%write out images
	bids(subj).func(1).results(11).name='cvr magnitude/phase maps from sinCVR - 5 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(11).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-ev-mag'];
	save_avw(sincvr4min_mag,bids(subj).func(1).results(11).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(11).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-ev-magvar'];
	save_avw(sincvr4min_magvar,bids(subj).func(1).results(11).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(11).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr4min-ev-magrsd'];
	save_avw(sincvr4min_magrsd,bids(subj).func(1).results(11).cvr_magrsd,'f',scales');

	%using 3 mins data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(12).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(12).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(12).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(12).feat 'stats/sigmasquareds']);
	
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

	bids(subj).func(1).results(12).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(12).petco2basevar=sinstd(3);
	bids(subj).func(1).results(12).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(12).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr3min_mag=cope./mean_func;
	sincvr3min_mag(isnan(sincvr3min_mag))=0;
	
	sincvr3min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr3min_magvar(isnan(sincvr3min_magvar))=0;
	
	sincvr3min_magrsd=sqrt(sincvr3min_magvar)./sincvr3min_mag;
	
	%write out images 
	bids(subj).func(1).results(12).name='cvr magnitude/phase maps from sinCVR - 3 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(12).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-ev-mag'];
	save_avw(sincvr3min_mag,bids(subj).func(1).results(12).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(12).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-ev-magvar'];
	save_avw(sincvr3min_magvar,bids(subj).func(1).results(12).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(12).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr3min-ev-magrsd'];
	save_avw(sincvr3min_magrsd,bids(subj).func(1).results(12).cvr_magrsd,'f',scales');
		
	%using 2 mins data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(13).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(13).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(13).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(13).feat 'stats/sigmasquareds']);
	
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

	bids(subj).func(1).results(13).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(13).petco2basevar=sinstd(3);
	bids(subj).func(1).results(13).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(13).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr2min_mag=cope./mean_func;
	sincvr2min_mag(isnan(sincvr2min_mag))=0;
	
	sincvr2min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr2min_magvar(isnan(sincvr2min_magvar))=0;
	
	sincvr2min_magrsd=sqrt(sincvr2min_magvar)./sincvr2min_mag;
	
	%write out images 
	bids(subj).func(1).results(13).name='cvr magnitude/phase maps from sinCVR - 2 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(13).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-ev-mag'];
	save_avw(sincvr2min_mag,bids(subj).func(1).results(13).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(13).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-ev-magvar'];
	save_avw(sincvr2min_magvar,bids(subj).func(1).results(13).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(13).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr2min-ev-magrsd'];
	save_avw(sincvr2min_magrsd,bids(subj).func(1).results(13).cvr_magrsd,'f',scales');
		
	%using 1 min data
	
	%read in source images
	[cope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(14).feat 'stats/cope1']);
	[mean_func dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(14).feat 'mean_func']);
	[varcope dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(14).feat 'stats/varcope1']);
	[sigmasqd dims scales bpp endian]=read_avw([bids(subj).func(1).analysis(14).feat 'stats/sigmasquareds']);
	
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

	bids(subj).func(1).results(14).petco2base=sina(3)-sindeltapetco2./2;
	bids(subj).func(1).results(14).petco2basevar=sinstd(3);
	bids(subj).func(1).results(14).petco2delta=sindeltapetco2;
	bids(subj).func(1).results(14).petco2deltavar=sinvardeltapetco2;
		
	%estimate magnitude and errors
	sincvr1min_mag=cope./mean_func;
	sincvr1min_mag(isnan(sincvr1min_mag))=0;
	
	sincvr1min_magvar=(cope./mean_func).^2.*(varcope./cope.^2+sigmasqd./mean_func.^2);
	sincvr1min_magvar(isnan(sincvr1min_magvar))=0;
	
	sincvr1min_magrsd=sqrt(sincvr1min_magvar)./sincvr1min_mag;
	
	%write out images 
	bids(subj).func(1).results(14).name='cvr magnitude/phase maps from sinCVR - 1 mins data';
	s=regexp(bids(subj).func(1).fname,'/');
	bids(subj).func(1).results(14).cvr_mag=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-ev-mag'];
	save_avw(sincvr1min_mag,bids(subj).func(1).results(14).cvr_mag,'f',scales');
	
	bids(subj).func(1).results(14).cvr_magvar=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-ev-magvar'];
	save_avw(sincvr1min_magvar,bids(subj).func(1).results(14).cvr_magvar,'f',scales');
	
	bids(subj).func(1).results(14).cvr_magrsd=[bids(subj).dir 'derivatives/' bids(subj).name '/results/' bids(subj).func(1).fname(s(end)+1:end) '_sincvr1min-ev-magrsd'];
	save_avw(sincvr1min_magrsd,bids(subj).func(1).results(14).cvr_magrsd,'f',scales');
