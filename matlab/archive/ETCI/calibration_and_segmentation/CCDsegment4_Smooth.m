function imageFull = CCDsegment4_Smooth(opts, imageOriginal)
%function imageFull = CCDsegment4_Smooth(opts, imageOriginal)
%
% Operated by CCDsegment4.

imageFull = filter2(opts.smoothingKernel,imageOriginal,'same');