
myVideo = VideoWriter('surfplots'); %open video file
myVideo.FrameRate = 2;  %can adjust this, 5 - 10 works well for me
open(myVideo)
p=[];
X=1:200;
Y=1:200;
[X,Y]=meshgrid(X,Y);
%% Plot in a loop and grab frames
for i=1:1:20
    for j=1:6
        A=data(j);
        subplot(2,3,j);
        surf(X,Y,A.hk(:,:,i),'LineStyle','none')
        zlim([-2 2])
        
       
        %biglim=get(p(j),'clim');
        %set(p(j),'clim',biglim);
        title(['t = ',num2str(i)])
        set(gcf, 'Position', get(0, 'Screensize'))
        colorbar
        caxis([-1 2])
    end
    
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
end
close(myVideo)
