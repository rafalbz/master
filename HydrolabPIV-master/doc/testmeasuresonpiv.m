if(~exist('piv6','var'))
    im1 = imread('imA.png');
    im2 = imread('imB.png');

    % Phase-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@phasecorr,true,@(x) x,'savepeaks',true); 
    piv1 = normalpass([],im1,[],im2,[],opt);

    % mqd
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedmqd,true,@(x) -x,'savepeaks',true); 
    piv2 = normalpass([],im1,[],im2,[],opt);

    % Cross-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedcc,true,@(x) x,'savepeaks',true); 
    piv3 = normalpass([],im1,[],im2,[],opt);

    % Normalized mad
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmdjold,1,@(x) -x,'savepeaks',true); 
    piv4b = normalpass([],im1,[],im2,[],opt);

    % Normalized mqd
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskednmqd,true,@(x) -x,'savepeaks',true); 
    piv5 = normalpass([],im1,[],im2,[],opt);

    % Normalized cross-corr
    opt = setpivopt('range',[-8 8 -8 8],'subwindow',24,24,.50,'measure',@maskedncc,true,@(x) x,'savepeaks',true); 
    piv6 = normalpass([],im1,[],im2,[],opt);
end

figure
cx = [-2.5 3.5];
subplot(2,3,1)
imagesc(piv1.U)
axis off
caxis(cx);
title('PC');

subplot(2,3,2)
imagesc(piv2.U)
axis off
caxis(cx);
title('MQD');

subplot(2,3,3)
imagesc(piv3.U)
axis off
caxis(cx);
title('CC');

subplot(2,3,4)
imagesc(piv4.U)
axis off
caxis(cx);
title('NMAD');

subplot(2,3,5)
imagesc(piv5.U)
caxis(cx);
title('NMQD');
axis off

subplot(2,3,6)
imagesc(piv6.U)
axis off
caxis(cx);
title('NCC');
print -depsc2 diffmeasures_velocity.eps

figure
subplot(2,3,1)
imagesc(piv1.F(:,:,end/2,end/2))
axis off
%caxis(cx);
title('PC');

subplot(2,3,2)
imagesc(-piv2.F(:,:,end/2,end/2))
axis off
%caxis(cx);
title('MQD');

subplot(2,3,3)
imagesc(piv3.F(:,:,end/2,end/2))
axis off
%caxis(cx);
title('CC');

subplot(2,3,4)
imagesc(-piv4.F(:,:,end/2,end/2))
axis off
caxis([-2 2]);
title('NMAD');

subplot(2,3,5)
imagesc(-piv5.F(:,:,end/2,end/2))
axis off
caxis([-2 2]);
title('NMQD');

subplot(2,3,6)
imagesc(piv6.F(:,:,end/2,end/2))
axis off
caxis([-1 1]);
title('NCC');
print -depsc2 diffmeasures_peaks.eps

