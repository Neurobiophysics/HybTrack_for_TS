function devmaxim = devmaxstack(PathFileName, frmmax)
im = imread(PathFileName,'tif',1);
im2 = zeros(size(im,1),size(im,2));
for i=1:frmmax
im = double(imread(PathFileName,'tif',i));
errdot = im>=10000;
avgim = mean(mean(im));
im(errdot) = avgim;
devim = im-avgim;
dsrim = (devim.^2).^(0.5);
im2 = cat(3,im2,dsrim);
end
im2= max(im2,[],3);
devmaxim=im2;