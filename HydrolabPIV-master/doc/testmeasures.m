if(~exist('piv4b','var'))
    im1 = imread('imA.png');
    im2 = imread('imB.png');
    [x,y] = meshgrid(-255.5:255.5);
    mask = (x.^2+y.^2)<256.^2 | x < 0 | y < 0;
    im1(~mask) = 0;
    im2(~mask) = 0;
    
    

    H = 32;
    W = 32;
    O = .75;
    % Phase-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@phasecorr,true,@(x) x,'savepeaks',true); 
    piv1a = normalpass([],im1,mask,im2,mask,opt);
    piv1b = normalpass(piv1a,im1,mask,im2,mask,opt);

    % Cross-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@standardcc,true,@(x) x,'savepeaks',true); 
    piv2a = normalpass([],im1,mask,im2,mask,opt);
    piv2b = normalpass(piv2a,im1,mask,im2,mask,opt);
    
    % Masked Cross-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@maskedcc,true,@(x) x,'savepeaks',true); 
    piv3a = normalpass([],im1,mask,im2,mask,opt);
    piv3b = normalpass(piv3a,im1,mask,im2,mask,opt);
    
    % Normalized cross-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@maskedncc,true,@(x) x,'savepeaks',true); 
    piv4a = normalpass([],im1,mask,im2,mask,opt);
    piv4b = normalpass(piv4a,im1,mask,im2,mask,opt);
            
%     % masked mqd
%     opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@maskedmqd,true,@(x) -x,'savepeaks',true); 
%     piv5 = normalpass([],im1,mask,im2,mask,opt);
%     
%     
%     % Normalized mad
%     opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@maskednmdj,1,@(x) -x,'savepeaks',true); 
%     piv6 = normalpass([],im1,mask,im2,mask,opt);
% 
%     % Normalized mqd
%     opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@maskednmqd,true,@(x) -x,'savepeaks',true); 
%     piv7 = normalpass([],im1,mask,im2,mask,opt);
% 
%     % Normalized cross-corr
%     opt = setpivopt('range',[-8 8 -8 8],'subwindow',H,W,O,'measure',@maskedncc,true,@(x) x,'savepeaks',true); 
%     piv8 = normalpass([],im1,mask,im2,mask,opt);
end
if(~exist('exactU','var'))
    load ~/face_460844/syntheticpiv/sizetest2/test1B/matpiv2/supermu_mp2.mat exactU
    U = exactU(:,:,1,1);
    V = exactU(:,:,1,2);
end
idx = 1:2:65;

figure;
cx = [-2.5 3.5];
subplot(2,3,1)
[U2,V2] = replaceoutliers(piv2a,mask);
%imagesc(sqrt((U2-U).^2 + (V2-V).^2));
imagesc(U2);
axis off
caxis(cx);
title('CC');
%colorbar

subplot(2,3,2)
[U3,V3] = replaceoutliers(piv3a,mask);
%imagesc(sqrt((U3-U).^2 + (V3-V).^2));
imagesc(U3);
axis off
caxis(cx);
title('Masked CC');
%colorbar

subplot(2,3,3)
[U4,V4] = replaceoutliers(piv4a,mask);
imagesc(U4);
axis off
caxis(cx);
title('Masked NCC');
%colorbar

cx = [-1 1]*.5;
subplot(2,3,4)
[U2b,V2b] = replaceoutliers(piv2b,mask);
%imagesc(sqrt((U2-U).^2 + (V2-V).^2));
imagesc(U2-U);
axis off
caxis(cx);
title('CC');
%colorbar

subplot(2,3,5)
[U3b,V3b] = replaceoutliers(piv3b,mask);
%imagesc(sqrt((U3-U).^2 + (V3-V).^2));
imagesc(U3-U);
axis off
caxis(cx);
title('Masked CC');
%colorbar

subplot(2,3,6)
[U4b,V4b] = replaceoutliers(piv4b,mask);
imagesc(sqrt((U4-U).^2 + (V4-V).^2));
axis off
caxis(cx);
title('Masked NCC');
%colorbar
% 
% figure
% cx = [-2.5 3.5];
% subplot(2,3,1)
% imagesc(piv1.U)
% axis off
% caxis(cx);
% title('PC');
% 
% subplot(2,3,2)
% imagesc(piv2.U)
% axis off
% caxis(cx);
% title('MQD');
% 
% subplot(2,3,3)
% imagesc(piv3.U)
% axis off
% caxis(cx);
% title('CC');
% 
% subplot(2,3,4)
% imagesc(piv4.U)
% axis off
% caxis(cx);
% title('NMAD');
% 
% subplot(2,3,5)
% imagesc(piv5.U)
% caxis(cx);
% title('NMQD');
% axis off
% 
% subplot(2,3,6)
% imagesc(piv6.U)
% axis off
% caxis(cx);
% title('NCC');
% print -depsc2 diffmeasures_velocity.eps
% 
% figure
% subplot(2,3,1)
% imagesc(piv1.F(:,:,end/2,end/2))
% axis off
% %caxis(cx);
% title('PC');
% 
% subplot(2,3,2)
% imagesc(-piv2.F(:,:,end/2,end/2))
% axis off
% %caxis(cx);
% title('MQD');
% 
% subplot(2,3,3)
% imagesc(piv3.F(:,:,end/2,end/2))
% axis off
% %caxis(cx);
% title('CC');
% 
% subplot(2,3,4)
% imagesc(-piv4.F(:,:,end/2,end/2))
% axis off
% caxis([-2 2]);
% title('NMAD');
% 
% subplot(2,3,5)
% imagesc(-piv5.F(:,:,end/2,end/2))
% axis off
% caxis([-2 2]);
% title('NMQD');
% 
% subplot(2,3,6)
% imagesc(piv6.F(:,:,end/2,end/2))
% axis off
% caxis([-1 1]);
% title('NCC');
% print -depsc2 diffmeasures_peaks.eps

