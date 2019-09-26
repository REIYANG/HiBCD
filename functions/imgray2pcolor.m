function rgb = imgray2pcolor(gim, map, n)
% IMGRAY2PSEUDOCOLOR transform a gray image to pseudocolor image
%   GIM is the input gray image data
%   MAP is the colormap already defined in MATLAB, for example:
%      'Jet','HSV','Hot','Cool','Spring','Summer','Autumn','Winter','Gray',
%      'Bone','Copper','Pink','Lines'
%   N specifies the size of the colormap 
%   rgb is the output COLOR image data
%
% Main codes stolen from:
%       http://www.alecjacobson.com/weblog/?p=1655
%       %% rgb = ind2rgb(gray2ind(im,255),jet(255));                      %
%                                                                           


[nr,nc,nz] = size(gim);
rgb = zeros(nr,nc,3);

if ( ~IsValidColormap(map) )
    disp('Error in ImGray2Pseudocolor: unknown colormap!');
elseif (~(round(n) == n) || (n < 0))
    disp('Error in ImGray2Pseudocolor: non-integer or non-positive colormap size');
else
    fh = str2func(ExactMapName(map));
    rgb = ind2rgb(gray2ind(gim,n),fh(n));
    rgb = uint8(rgb*255);
end

if (nz == 3)
    rgb = gim;
    disp('Input image has 3 color channel, the original data returns');
end

end

function y = IsValidColormap(map)

y = strncmpi(map,'Jet',length(map)) | strncmpi(map,'HSV',length(map)) |...
    strncmpi(map,'Hot',length(map)) | strncmpi(map,'Cool',length(map)) |...
    strncmpi(map,'Spring',length(map)) | strncmpi(map,'Summer',length(map)) |...
    strncmpi(map,'Autumn',length(map)) | strncmpi(map,'Winter',length(map)) |...
    strncmpi(map,'Gray',length(map)) | strncmpi(map,'Bone',length(map)) |...
    strncmpi(map,'Copper',length(map)) | strncmpi(map,'Pink',length(map)) |...
    strncmpi(map,'Lines',length(map));
end

function emapname = ExactMapName(map)

if strncmpi(map,'Jet',length(map))
    %emapname = 'Jet';
    emapname = 'jet';
elseif strncmpi(map,'HSV',length(map))
    %emapname = 'HSV';
    emapname = 'hsv';
elseif strncmpi(map,'Hot',length(map))
    %emapname = 'Hot';
    emapname = 'hot';
elseif strncmpi(map,'Cool',length(map))
    %emapname = 'Cool';
    emapname = 'cool';
elseif strncmpi(map,'Spring',length(map))
    %emapname = 'Spring';
    emapname = 'spring';
elseif strncmpi(map,'Summer',length(map))
    %emapname = 'Summer';
    emapname = 'summer';
elseif strncmpi(map,'Autumn',length(map))
    %emapname = 'Autumn';
    emapname = 'autumn';
elseif strncmpi(map,'Winter',length(map))
    %emapname = 'Winter';
    emapname = 'winter';
elseif strncmpi(map,'Gray',length(map))
    %emapname = 'Gray';
    emapname = 'gray';
elseif strncmpi(map,'Bone',length(map))
    %emapname = 'Bone';
    emapname = 'bone';
elseif strncmpi(map,'Copper',length(map))
    %emapname = 'Copper';
    emapname = 'copper';
elseif strncmpi(map,'Pink',length(map))
    %emapname = 'Pink';
    emapname = 'pink';
elseif strncmpi(map,'Lines',length(map))
    %emapname = 'Lines';
    emapname = 'lines';
end

end 