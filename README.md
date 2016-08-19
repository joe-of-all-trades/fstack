# fstack

FSTACK merging images of mutiple focal planes into one in-focus image. This script is a polished version of the codes posted at this link : http://www.magiclantern.fm/forum/index.php?topic=11886.0 It works on the sample dataset found at the link as well. For each dataset, parameter has to be tuned carefully. 

Input must be a cell array of images taken at different focal planes, with each cell containing a single color or grayscale image.

edofimg = fstack(img) merges img, an img array containing grayscale or
color images acquired at mutiple focal distance, into one all-in-focus
image.

edofimg = fstack(img, option, value) performs image merging using user 
specified option value pair. 
e.g. edofimg = fstack(img,'logSize',15);

[edofimg, fmap] = fstack(img) returns an all-in-focus image and a image
of focal planes from which a particular pixel is extracted.

[edofimg, fmap, logrespone] = fstack(img) returns an all-in-focus
imgage, image of focal planes, and maximum logrespone image. The
logresponse image can be used to set a logresponse threshold to improve
result.

available options :

'logsize'    : size of the LoG (laplacian of gaussian) filter used for 
            detecting pixel in focus, default is 13

'logstd'     : standard deviation for LoG filter, default is 2
'dilatesize' : size of structure element used to smooth the focus
            detection result, default is 31
'blendsize'  : size of the Guassian filter for bleding pixels taken from
            different focal planes, default is 31
'blendstd'   : standard deviation of the Gaussian filter for blending
            pixels from different planes, default is 5
'logthreshold' : threshold for logresponse, default is 0

Version 1.0
Copyright: Chao-Yuan Yeh, 2016
