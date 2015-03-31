function callTrackTreatP2(N)

str=['load("tmp',num2str(N),'.mat");'];
str=rmquote(str);
eval(str);

eval(['[TTa,TTp,TTe,TThw]=TrackTreatP3a(out',num2str(N),',3,3);']);
str=['save("TrackerOuts',num2str(N),'_3.mat","TTa","TTp","TTe","TThw");'];
str=rmquote(str);
eval(str);
