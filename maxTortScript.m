

maxTort = 0;
maxiV = 0;
maxiB = 0;

minTort = Inf;
miniV = 0;
miniB = 0;

for iV = 1:numel(V)
    
   for iB = 1:numel(V{iV}.Branching.Branches)
       
       Tort = V{iV}.Branching.Branches{iB}.Tortuosity;
       
       if (Tort > maxTort && size(V{iV}.Branching.Branches{iB}.Points,1) > 30)
           maxTort = Tort;
           maxiV = iV;
           maxiB = iB;
           maxPoints = V{iV}.Branching.Branches{iB}.SmoothedPoints;
       end
       
       if (Tort < minTort && size(V{iV}.Branching.Branches{iB}.Points,1) > 30)
           minTort = Tort;
           miniV = iV;
           miniB = iB;
           minPoints = V{iV}.Branching.Branches{iB}.SmoothedPoints;
       end
       
   end
    
   
end
bMovie = 0;
figure(1)
PlotBranch(V,maxiV,maxiB);
axis equal;
title('Max Tortuosity');
name = 'MaxTortMovie';

 if(bMovie)
        for iA = 1:360
            view(iA,15);
            pause(0.01);
                 if iA == 1;
                          writerObj = VideoWriter([name,'.avi']);
                          open(writerObj);
                          frame = getframe(1);
                          writeVideo(writerObj,frame);
                  else
                          frame = getframe(1);
                          writeVideo(writerObj,frame);
                 end
        end
end

figure(2)
PlotBranch(V,miniV,miniB);
axis equal;
title('Min Tortuosity');
name = 'MinTortMovie';

if(bMovie)
    for iA = 1:360
        view(iA,15);
        pause(0.01);
             if iA == 1;
                      writerObj = VideoWriter([name,'.avi']);
                      open(writerObj);
                      frame = getframe(2);
                      writeVideo(writerObj,frame);
              else
                      frame = getframe(2);
                      writeVideo(writerObj,frame);
             end
    end
end

close(writerObj);