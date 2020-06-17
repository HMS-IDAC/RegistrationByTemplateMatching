I = im2double(imread('cell.tif'));
J = imtranslate(imrotate(I,5,'crop'),[5 10]);
tform = rtm(J,I,'regType','rigid','angleSet',-15:5:15);
TJ = imwarp(J,tform,'OutputView',imref2d(size(I)));
subplot(1,3,1), imshow(I), title('fixed')
subplot(1,3,2), imshow(J), title('moving')
subplot(1,3,3), imshow(TJ), title('registered')