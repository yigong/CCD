function callTrackTreatP2(N)

str=['load("tmp',num2str(N),'.mat");'];
str=rmquote(str);
eval(str);
eval(['[NPL]=get_locs_and_starts_2014(out',num2str(N),');']);
NT=length(NPL)
eval(['[TTa,TTp,TTe,TThw]=TrackTreatP2(out',num2str(N),');']);
str=['save("TrackerOuts',num2str(N),'.mat","NT","TTa","TTp","TTe","TThw");'];
str=rmquote(str);
eval(str);


