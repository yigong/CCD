function callTrackTreatP3(N,a)
% N: file name suffix
% a: True alpha
str=['load("/Users/Yigong/Research/photon_diagnostics/data/sampleOnCCD/tmp',num2str(N),'.mat");'];
str=rmquote(str);
eval(str);
str=['[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(out',num2str(N),');']
eval(['[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(out',num2str(N),');']);
NT=length(no_primary_locs)

['start 1st in tmp', num2str(N),'.mat']
eval(['[alpha,alphaTrue,energyMeasured,wrongEnd,beta,ridgeWidth]=TrackTreatP3R(out',num2str(N),',a,3,1);']);
str=['save("TrackerOuts',num2str(N),'_1.mat","NT","alpha","alphaTrue","beta","energyMeasured","wrongEnd","ridgeWidth");'];
str=rmquote(str);
eval(str);

['start 2nd in tmp', num2str(N), '.mat']
eval(['[alpha,alphaTrue,energyMeasured,wrongEnd,beta,ridgeWidth]=TrackTreatP3R(out',num2str(N),',a,3,2);']);
str=['save("TrackerOuts',num2str(N),'_2.mat","NT","alpha","alphaTrue","beta","energyMeasured","wrongEnd","ridgeWidth");'];
str=rmquote(str);
eval(str);

['start 3rd in tmp', num2str(N),'.mat']
eval(['[alpha,alphaTrue,energyMeasured,wrongEnd,beta,ridgeWidth]=TrackTreatP3R(out',num2str(N),',a,3,3);']);
str=['save("TrackerOuts',num2str(N),'_3.mat","NT","alpha","alphaTrue","beta","energyMeasured","wrongEnd","ridgeWidth");'];
str=rmquote(str);
eval(str);
