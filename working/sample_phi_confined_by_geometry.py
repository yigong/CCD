import numpy as np

windowX_sampled = np.array([0]) # in mm
windowZ_sampled = np.array([0]) # in mm

rowPixel = 1500.
colPixel = 1500.

yEdge = 2.97
zEdge = 1.0

pxSize = 10.5/1000 # in mm
CCDlowerLeft = {"x": -1753.5*pxSize, "y": yEdge+0.5*pxSize, "z": 0} # in mm

xPixel = CCDlowerLeft["x"] + colPixel * pxSize
yPixel = CCDlowerLeft["y"] + rowPixel * pxSize
zPixel = 0
alphaRadius = np.abs(np.arctan2(yPixel, xPixel) - np.pi/2)
cosAlpha = np.cos(alphaRadius)

windowLower = {"x": 0, "y": 0, "z": zEdge} # in mm
xWindow = windowX_sampled + windowLower["x"]
yWindow = np.zeros_like(xWindow)
zWindow = windowZ_sampled + windowLower["z"]
cosBeta = np.cos(np.arctan2(zWindow, np.sqrt(xPixel**2 + yPixel**2)))

cosPhi = yPixel/np.sqrt((xPixel-xWindow)**2 + \
                        (yPixel-yWindow)**2 + (zPixel-zWindow)**2)
cosPhi_check = cosAlpha * cosBeta # check the cosine relation

mc2 = 511.








