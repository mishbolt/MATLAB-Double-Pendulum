
myVideo = VideoWriter('cplots'); %open video file
myVideo.FrameRate = 2;  %can adjust this, 5 - 10 works well for me
open(myVideo)
p=[];
%% Plot in a loop and grab frames
for i=1:1:20
    for j=1:6
        A=data(j);
        subplot(2,3,j);
        contourf(A.hk(:,:,i),'LineColor','none')
        caxis([-1 2])
        %biglim=get(p(j),'clim');
        %set(p(j),'clim',biglim);
        title(['t = ',num2str(i)])
        set(gcf, 'Position', get(0, 'Screensize'))
        colorbar
    end
    
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
end
close(myVideo)
