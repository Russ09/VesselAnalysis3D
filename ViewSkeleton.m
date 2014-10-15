function [] = ViewSkeleton( Skeleton, Settings )
%Visualise Skeleton of vessels
if iscell(Skeleton)
    
    figure(1)
    VesselsToDo = 1:numel(Skeleton);
    cols = colormap;
    
    if nargin > 1
        if isfield(Settings,'Vessels')
            VesselsToDo = Settings.Vessels;
        end
        
        if isfield(Settings,'Colormap')
            Colormap = Settings.Colormap
        end
        
    end
    
    for iV = VesselsToDo
        
        CurrVes = Skeleton{iV};
        nBranches = numel(CurrVes.Branching.Branches);
        cols = colormap(jet(nBranches));
        
        for iB = 1:nBranches
            
            Branch = CurrVes.Branching.Branches{iB};
            Points = Branch.SmoothedPoints;
            
            plot3(Points(:,1),Points(:,2),Points(:,3),'color',cols(iB,:),'LineWidth',2);
            hold on
        end
    end
    
    hold off
    axis equal
    bMovie = 1;
    name = 'SkeletonMovie';
    
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
        
    
else
    radius = 3;
    stencil = ones(radius,radius,radius);

    DilateSkel = imdilate(Skeleton,stencil);

    isosurface(DilateSkel)
    axis equal
    lighting gouraud
    
end

end

