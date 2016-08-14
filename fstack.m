function [edofimg, varargout] = fstack(img,varargin)
%FSTACK merging images of mutiple focal planes into one in-focus image.
%
%  edofimg = fstack(img) merges img, an img array containing grayscale or
%  color images acquired at mutiple focal distance, into one all-in-focus
%  image.
%  
%  edofimg = fstack(img, option, value) performs image merging using user 
%  specified option value pair. 
%  e.g. edofimg = fstack(img,'logSize',15);
%   
%  [edofimg, fmap] = fstack(img) returns an all-in-focus image and a image
%  of focal planes from which a particular pixel is extracted.
%
%  [edofimg, fmap, logrespone] = fstack(img) returns an all-in-focus
%  imgage, image of focal planes, and maximum logrespone image. The
%  logresponse image can be used to set a logresponse threshold to improve
%  result.
%  
%  available options :
%
%  'logsize'    : size of the LoG (laplacian of gaussian) filter used for 
%                 detecting pixel in focus, default is 13
%
%  'logstd'     : standard deviation for LoG filter, default is 2
%  'dilatesize' : size of structure element used to smooth the focus
%                 detection result, default is 31
%  'blendsize'  : size of the Guassian filter for bleding pixels taken from
%                 different focal planes, default is 31
%  'blendstd'   : standard deviation of the Gaussian filter for blending
%                 pixels from different planes, default is 5
%  'logthreshold' : threshold for logresponse, default is 0
%  
%  Version 1.0
%  Copyright: Chao-Yuan Yeh, 2016
%  
% Cell array allows more flexible indexing. Some of the logic indexing
% method used here won't work with multi-dimensional array. Logical
% indexing dramatically improves execution speed. 
if ~iscell(img)
    error('Input needs to be cell array of images.')
end

filtersize = checkparam('logsize', 13);
% Keep in mind that the size of LoG filter should ideally be the smallest
% odd number that is larger than 6 standard deviation so filter effect
% won't be truncated. Larger filter size is unnecessary. 
logstd = checkparam('logstd', 2);
dilatesize = checkparam('dilatesize', 31);
blendsize = checkparam('blendsize', 31);
blendstd = checkparam('blendstd', 5);
logthreshold = checkparam('logthreshold', 0);

function val = checkparam(param, defaultval)
    global varargin
    if any(strcmpi(varargin, param))
        val = varargin{find(strcmpi(varargin, param))+1};
    else
        val = defaultval;
    end
end

if (numel(size(img{1}))>2)
    imgcopy = img;
    colorimg = true;
    for ii = 1:length(imgcopy)
        img{ii} = rgb2gray(img{ii});
    end
else
    colorimg = false;
end

% make brightness of all images equal
if ~colorimg
    avg1 = mean2(img{1});
    for ii = 2 : length(img)
        avgcur = mean2(img{ii});
        img{ii} = img{ii} + avg1 - avgcur;
    end
else
    avg1 = mean(imgcopy{1}(:));
    for ii = 2:length(imgcopy)
    avgcur = mean(imgcopy{ii}(:));
    imgcopy{ii} = imgcopy{ii} + ceil(avg1 - avgcur);  
    end
end

imgfiltered = cell(size(img));

logfilter = fspecial('log', [filtersize filtersize], logstd);
se = strel('ball', dilatesize, dilatesize);

for ii = 1:length(img)    
    imgfiltered{ii} = imfilter(single(img{ii}), logfilter);
    % Note that LoG detects border with zero-crossing. It might be worth
    % playing with the following. 
    % imgfiltered{ii} = -imfilter(single(img{ii}), logfilter);
    % imgFiltered{ii} = abs(imgFiltered{ii});
    imgfiltered{ii} = imdilate(imgfiltered{ii}, se, 'same');
end

fmap = ones(size(img{1}), 'single');
logresponse = zeros(size(img{1}), 'single') + logthreshold;

% Look for focal plane that has the largest LoG response (pixel in focus). 
for ii = 1:length(img)
    index = imgfiltered{ii} > logresponse;
    logresponse(index) = imgfiltered{ii}(index);
    fmap(index) = ii;
end

% Smooth focal plane image
fmap = imfilter(fmap,fspecial('gaussian', [blendsize blendsize], blendstd));
fmap(fmap < 1) = 1;

if ~colorimg
    edofimg = img{1};
else
    edofimg = imgcopy{1};
end

% Extract in-focus pixel from every image. 
for ii = 1:length(img)
    index = fmap == ii;
    if ~colorimg
        edofimg(index) = img{ii}(index);
    else
        index = repmat(index,[1 1 3]);
        edofimg(index) = imgcopy{ii}(index);
    end
end

% Blend different focal planes
for ii = 1:length(img)-1
    index = fmap > ii & fmap < ii+1;
    if ~colorimg
        edofimg(index) = (fmap(index) - ii).*single(img{ii+1}(index)) + ...
            (ii+1-fmap(index)).*single(img{ii}(index));
    else
        index = repmat(index,[1 1 3]);
        fmap_c = repmat(fmap, [1 1 3]);
        edofimg(index) = ( fmap_c(index) - ii).*single(imgcopy{ii+1}(index)) + ...
            (ii + 1 - fmap_c(index)).*single(imgcopy{ii}(index));
    end
end

if nargout >= 2
    varargout{1} = fmap;
end

if nargout == 3
    varargout{2} = logresponse;
end

end


    