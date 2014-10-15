function [] = SpinFrameMovie(name,el)

bAvi = 1;
bGif = 0;

if(nargin < 2)
    el = 15;
end

if(bAvi)
    for iA = 1:360
                view(90,iA)
                hFig = figure(1);
                set(hFig,'Position',[50,50,1000,1000]);
                pause(0.01)
                     if iA == 1
                              writerObj = VideoWriter([name,'.avi']);
                              open(writerObj);
                              frame = getframe(gcf);
                              writeVideo(writerObj,frame);
                      else
                              frame = getframe(gcf);
                              writeVideo(writerObj,frame);
                     end
            end

end

if (bGif)
    
       for iA = 1:360
                view(iA,el);
                
                pause(0.01);
                     if iA == 1
                              frame = getframe();
                              im = frame2im(frame);
                              [imind,cm] = rgb2ind(im,256);
                              imwrite(imind,cm,name,'gif', 'Loopcount',inf);
                      else
                              frame = getframe();
                              im = frame2im(frame);
                              [imind,cm] = rgb2ind(im,256);
                              imwrite(imind,cm,name,'gif','WriteMode','append');
                     end
       end

end 
    
    