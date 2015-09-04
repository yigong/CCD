basedirin = '/Users/Yigong/Research/photon_diagnostics/data/08_08_14_background/cropped/';
dirlist{200, 1} = []; 
flist = dir([basedirin, '*.fit'])
filesuffix = '_median_subtracted'
basedirout = '/Users/Yigong/Research/photon_diagnostics/data/08_08_14_background/median_subtracted/';
gain = ones(1,4);
replaceflag = false

CCDcal_median2(basedirin,dirlist,flist,filesuffix,basedirout,gain,replaceflag)