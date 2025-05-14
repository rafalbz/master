

if(~exist('piv5n','var'))
    % loss of particles test?
    %im1 = im2double(imread('im50x000001A.png'));
    %im2 = im2double(imread('im25x000001B.png'));
    im1 = im2double(imread('imA.png'));
    im2 = im2double(imread('imB.png'));
    
    %noise = 1/4; 
    noise = 1/4096;
    im1 = im1 + noise*rand(512);
    im2 = im2 + noise*rand(512);    
    
    
    q = 8;
    %im1 = round(q*im1)/q + noise*rand(512);
    %im2 = round(q*im2)/q + noise*rand(512);
    idx = (-15.5:15.5) + 250.5;
    A = im1(idx,idx);
    B = im2(idx,idx);
    M = ones(size(A)); 
    
    % alpha = 0
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedmdj,0,@(x) -x,'savepeaks',true); 
    piv1 = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmdj,0,@(x) -x,'savepeaks',true); 
    piv1n = normalpass([],im1,[],im2,[],opt);
    F1 = maskedmdj(A,M,B,M,[],[],0);
    F1n = maskednmdj(A,M,B,M,[],[],0);
    
    % alpha = 1
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedmdj,1,@(x) -x,'savepeaks',true); 
    piv2 = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmdj,1,@(x) -x,'savepeaks',true); 
    piv2n = normalpass([],im1,[],im2,[],opt);
    F2 = maskedmdj(A,M,B,M,[],[],1);
    F2n = maskednmdj(A,M,B,M,[],[],1);
        
     % alpha = 2
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedmdj,2,@(x) -x,'savepeaks',true); 
    piv3 = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedmqd,true,@(x) -x,'savepeaks',true); 
    piv3b = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmdj,2,@(x) -x,'savepeaks',true); 
    piv3n = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmqd,true,@(x) -x,'savepeaks',true); 
    piv3nb = normalpass([],im1,[],im2,[],opt);
    
    F3 = maskedmdj(A,M,B,M,[],[],2);
    F3n = maskednmdj(A,M,B,M,[],[],2);
    F3b = maskedmqd(A,M,B,M);
    F3nb = maskednmqd(A,M,B,M);
    
    % alpha = Inf
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedmdj,inf,@(x) -x,'savepeaks',true); 
    piv4 = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmdj,inf,@(x) -x,'savepeaks',true); 
    piv4n = normalpass([],im1,[],im2,[],opt);
    F4 = maskedmdj(A,M,B,M,[],[],inf);
    F4n = maskednmdj(A,M,B,M,[],[],inf);
    
    % Cross-correlation
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedcc,true,@(x) x,'savepeaks',true); 
    piv5 = normalpass([],im1,[],im2,[],opt);
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedncc,true,@(x) x,'savepeaks',true); 
    piv5n = normalpass([],im1,[],im2,[],opt);
    F5 = maskedmdj(A,M,B,M,[],[],NaN);
    F5n = maskednmdj(A,M,B,M,[],[],NaN);
 
end

figure
cx = [-2.5 3.5];
subplot(4,5,1)
imagesc(piv1.U)
axis off equal
caxis(cx);
title('\alpha=0');

subplot(4,5,2)
imagesc(piv2.U)
axis off equal
caxis(cx);
title('\alpha=1');


subplot(4,5,3)
imagesc(piv3.U)
axis off  equal
caxis(cx);
title('\alpha=2');


subplot(4,5,4)
imagesc(piv4.U)
axis off equal
caxis(cx);
title('\alpha=\infty');


subplot(4,5,5)
imagesc(piv5.U)
axis off equal
caxis(cx);
title('CC');

subplot(4,5,6)
imagesc(piv1n.U)
axis off equal
caxis(cx);
title('\alpha=0');

subplot(4,5,7)
imagesc(piv2n.U)
axis off equal
caxis(cx);
title('\alpha=1');


subplot(4,5,8)
imagesc(piv3n.U)
axis off equal
caxis(cx);
title('\alpha=2');


subplot(4,5,9)
imagesc(piv4n.U)
axis off equal
caxis(cx);
title('\alpha=\infty');


subplot(4,5,10)
imagesc(piv5n.U)
axis off equal
caxis(cx);
title('NCC');

drawnow;
%print -depsc2 diffmeasures_velocity2.eps


idx = (-12:12) + 32; %+24;
figure
%cx = [-2.5 3.5];
subplot(4,5,11)
imagesc(F1(idx,idx))
axis off equal
%caxis(cx);
title('\alpha=0');

subplot(4,5,12)
imagesc(F2(idx,idx))
axis off equal
%caxis(cx);
title('\alpha=1');


subplot(4,5,13)
imagesc(F3(idx,idx))
axis off  equal
%caxis(cx);
title('\alpha=2');


subplot(4,5,14)
imagesc(F4(idx,idx))
axis off equal
%caxis(cx);
title('\alpha=\infty');


subplot(4,5,15)
imagesc(F5(idx,idx))
axis off equal
%caxis(cx);
title('CC');

cx = [0 2];
subplot(4,5,16)
imagesc(F1n(idx,idx))
axis off equal
%caxis(cx);
title('\alpha=0');

subplot(4,5,17)
imagesc(F2n(idx,idx))
axis off equal
%caxis(cx);
title('\alpha=1');


subplot(4,5,18)
imagesc(F3n(idx,idx))
axis off  equal
%caxis(cx);
title('\alpha=2');


subplot(4,5,19)
imagesc(F4n(idx,idx))
axis off equal
%caxis(cx);
title('\alpha=\infty');


subplot(4,5,20)
imagesc(F5n(idx,idx))
axis off equal
%caxis([-1 1]);
title('NCC');

drawnow;
%print -depsc2 diffmeasures_peaks2.eps

if(false)
scaleim = @(im) round(255*(im - min(im(:)))/(max(im(:))-min(im(:))) + 1);

imwrite(ind2rgb(scaleim(F1(idx,idx)),parula(256)),'measureimages/F1.png')
imwrite(ind2rgb(scaleim(F2(idx,idx)),parula(256)),'measureimages/F2.png')
imwrite(ind2rgb(scaleim(F3(idx,idx)),parula(256)),'measureimages/F3.png')
imwrite(ind2rgb(scaleim(F4(idx,idx)),parula(256)),'measureimages/F4.png')
imwrite(ind2rgb(scaleim(F5(idx,idx)),parula(256)),'measureimages/F5.png')
imwrite(ind2rgb(scaleim(F1n(idx,idx)),parula(256)),'measureimages/F1n.png')
imwrite(ind2rgb(scaleim(F2n(idx,idx)),parula(256)),'measureimages/F2n.png')
imwrite(ind2rgb(scaleim(F3n(idx,idx)),parula(256)),'measureimages/F3n.png')
imwrite(ind2rgb(scaleim(F4n(idx,idx)),parula(256)),'measureimages/F4n.png')
imwrite(ind2rgb(scaleim(F5n(idx,idx)),parula(256)),'measureimages/F5n.png')

scaleim = @(im) min(max(round(255*(im + 2.5)/(3.5+2.5) + 1),1),256);
imwrite(ind2rgb(scaleim(piv1.U),parula(256)),'measureimages/U1.png')
imwrite(ind2rgb(scaleim(piv2.U),parula(256)),'measureimages/U2.png')
imwrite(ind2rgb(scaleim(piv3.U),parula(256)),'measureimages/U3.png')
imwrite(ind2rgb(scaleim(piv4.U),parula(256)),'measureimages/U4.png')
imwrite(ind2rgb(scaleim(piv5.U),parula(256)),'measureimages/U5.png')
imwrite(ind2rgb(scaleim(piv1n.U),parula(256)),'measureimages/U1n.png')
imwrite(ind2rgb(scaleim(piv2n.U),parula(256)),'measureimages/U2n.png')
imwrite(ind2rgb(scaleim(piv3n.U),parula(256)),'measureimages/U3n.png')
imwrite(ind2rgb(scaleim(piv4n.U),parula(256)),'measureimages/U4n.png')
imwrite(ind2rgb(scaleim(piv5n.U),parula(256)),'measureimages/U5n.png')
end
