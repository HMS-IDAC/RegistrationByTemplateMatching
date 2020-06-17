function tform = rtm(J,I,varargin) % J: moving, I: fixed
% tform = rtm(J,I,varargin)
% registration by template matching
% registers moving image (J) w.r.t. fixed image (I)
% supports 'translation', 'rigid', and 'similarity' transforms
% see m-file for parameter details
% 
% ----- example (see RunMe.m) -----
% 
% I = im2double(imread('cell.tif'));
% J = imtranslate(imrotate(I,5,'crop'),[5 10]);
% tform = rtm(J,I,'regType','rigid','angleSet',-15:5:15);
% TJ = imwarp(J,tform,'OutputView',imref2d(size(I)));
% subplot(1,3,1), imshow(I), title('fixed')
% subplot(1,3,2), imshow(J), title('moving')
% subplot(1,3,3), imshow(TJ), title('registered')
%
% ----- reference -----
% The algorithm is described in Subsection 3.2 of
%   Marcelo Cicconet, David G. C. Hildebrand, and Hunter Elliott.
%   Finding Mirror Symmetry via Registration and Optimal Symmetric Pairwise Assignment of Curves.
%   IEEE ICCV, Detecting Symmetry in the Wild Workshop, 2017.
%   http://openaccess.thecvf.com/content_ICCV_2017_workshops/papers/w24/Cicconet_Finding_Mirror_Symmetry_ICCV_2017_paper.pdf
%
% Marcelo Cicconet, Dec 21 2017

ip = inputParser;
ip.addParameter('regType','similarity'); % options: 'translation', 'rigid', 'similarity'
ip.addParameter('scaleSet',0.75:0.125:1.25); % scales considered
ip.addParameter('angleSet',0:6:360-6); % rotations considered
ip.addParameter('maxDisp',Inf); % maximum distance of translations considered
ip.addParameter('boxSize',round(0.25*min(size(J)))); % template size for normxcorr2
ip.addParameter('nBoxSamples',100); % maximum number of template samples voting
ip.parse(varargin{:});
p = ip.Results;
regType = p.regType;
scaleSet = p.scaleSet;
angleSet = p.angleSet;
maxDisp = p.maxDisp;
boxSize = p.boxSize;
nBoxSamples = p.nBoxSamples;

if strcmp(regType,'similarity')
    allTforms = {};
    allSrmags = [];
    sclIndices = [];
    for iScale = 1:length(scaleSet)
        for i = 1:iScale
            fprintf('|');
        end
        for i = iScale+1:length(scaleSet)
            fprintf('.');
        end
        scale = scaleSet(iScale);
        scaleTform = affine2d([scale 0 0; 0 scale 0; 0 0 1]);
        RJ = imwarp(J,scaleTform);
        [tforms, srmags] = rrtm(RJ,I,angleSet,maxDisp,boxSize,nBoxSamples,10);
        allTforms = {allTforms{:} tforms{:}};
        allSrmags = [allSrmags srmags];
        sclIndices = [sclIndices iScale*ones(1,length(srmags))];
        for i = 1:length(scaleSet)
            fprintf('\b');
        end
    end
    [~,im] = max(allSrmags);
    scale = scaleSet(sclIndices(im));
    scaleTformT = [scale 0 0; 0 scale 0; 0 0 1];
    rigidTform = allTforms{im};
    tform = affine2d(rigidTform.T*scaleTformT);
elseif strcmp(regType,'rigid')
    tforms = rrtm(J,I,angleSet,maxDisp,boxSize,nBoxSamples,1);
    tform = tforms{1};
elseif strcmp(regType,'translation')
    tforms = rrtm(J,I,0,maxDisp,boxSize,nBoxSamples,1);
    tform = tforms{1};
else
    error('regType not recognized')
end

end

% -------------------------

function [tforms, srmags] = rrtm(J,I,angleSet,maxDisp,boxSize,nBoxSamples,maxNOutputs) % J: moving, I: fixed
% rigid registration by template matching

[numrows,numcols] = size(I);
rangles = angleSet;
rmags = zeros(1,length(rangles)); % likelihoods of most likely translation per angle
vs = zeros(length(rangles),2); % for every angle, the most likely translation

parfor iangle = 1:length(rangles)
    A = zeros(2*numrows,2*numcols);
    rotJ = imrotate(J,rangles(iangle),'crop');
    for index = 1:nBoxSamples
        w = boxSize;
        h = w;
        x0 = floor((size(rotJ,2)-w)*rand);
        y0 = floor((size(rotJ,1)-h)*rand);

        xcSubRotJ = x0+w/2;
        ycSubRotJ = y0+h/2;

        subRotJ = imcrop(rotJ,[x0 y0 w h]);
        if var(subRotJ(:)) > 0.001
            [ROI,~,mc] = locateSubset(subRotJ,I);
            if mc > 0.25
                xcI = ROI(1)+ROI(3)/2;
                ycI = ROI(2)+ROI(4)/2;

                v = [xcI ycI]-[xcSubRotJ ycSubRotJ]; % translation
                if norm(v) <= maxDisp
                    row = max(min(round(v(2)+numrows),2*numrows),1);
                    col = max(min(round(v(1)+numcols),2*numcols),1);
                    A(row,col) = A(row,col)+1;
                end
            end
        end
    end
    A = imgaussfilt(A,3);
    
    maxA = max(A(:));
    
    [r,c] = find(A == maxA);
    v(1) = c(1)-numcols;
    v(2) = r(1)-numrows;
    rmags(iangle) = maxA;
    vs(iangle,:) = v;
end

[srmags,iangles] = sort(rmags,'descend');
nTForms = min(length(rangles),maxNOutputs);
tforms = cell(1,nTForms);
for i = 1:nTForms
    iangle = iangles(i); % rotation
    
    v = vs(iangle,:); % translation

    arad = -rangles(iangle)/360*2*pi;

    % rotation with respect to center
    T1 = [eye(2) [-numcols/2; -numrows/2]; 0 0 1];
    T2 = [[cos(arad) -sin(arad); sin(arad) cos(arad)] [0; 0]; 0 0 1];
    T3 = [eye(2) [numcols/2; numrows/2]; 0 0 1];

    % translation
    T4 = [eye(2) v'; 0 0 1];

    % transform
    tform = affine2d((T4*T3*T2*T1)');
    
    tforms{i} = tform;
end

end

% -------------------------

function [ROI,c,mc] = locateSubset(ss,I)
% template matching

c = normxcorr2(ss,I);

mc = max(c(:));
[ypeak, xpeak] = find(c==mc);

yoffSet = ypeak(1)-size(ss,1);
xoffSet = xpeak(1)-size(ss,2);

ROI = [xoffSet+1, yoffSet+1, size(ss,2), size(ss,1)];

end