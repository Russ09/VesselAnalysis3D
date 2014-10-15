function [ output_args ] = PlotBranch( Vessels,nVessel,nBranch,bAll,options )


if(nargin >= 4 && bAll == 1)
    nVessel = 1:numel(Vessels)
end

if(nargin > 4 && isfield(options,'Colour'))
    Colour = options.Colour;
else
    Colour = [94,122,201]/255;
end

if(nargin > 4 && isfield(options,'bSkeleton'))
    bSkeleton = options.bSkeleton;
else
    bSkeleton = 0;
end

    
for iV = 1:length(nVessel)
    
    if(nargin >= 4 && bAll == 1)
        nBranch = 1:numel(Vessels{nVessel(iV)}.Branching.Branches);
    end
    
    if(isa(nBranch,'char')&&strcmp(nBranch,'all'))
        nBranch = 1:numel(Vessels{nVessel(iV)}.Branching.Branches);
    end
    
    for iB = 1:length(nBranch)
        Branch = Vessels{nVessel(iV)}.Branching.Branches{nBranch(iB)};

        Points = Branch.SmoothedPoints(1:end,:);
        
        if (isfield(Branch,'Up')&&size(Branch.Up,1)>0)
            Up = Branch.Up;
            Down = Branch.Down;
            Left = Branch.Left;
            Right = Branch.Right;

            if (isfield(Branch,'UpLeft')&&numel(Branch.UpLeft) == numel(Up))
                UpLeft = Branch.UpLeft;
                DownRight = Branch.DownRight;
                DownLeft = Branch.DownLeft;
                UpRight = Branch.UpRight;
                try
                    Boundary = cat(3,Up,UpLeft,Left,UpRight,Down,DownRight,Right,DownLeft,Up);
                catch
                    disp('Error');
                end
                    
            else
                Boundary = cat(3,Up,Down,Left,Right,Up);
            end

            plot3(Points(:,2),Points(:,1),Points(:,3),'o','MarkerEdgeColor','k','MarkerFaceColor',[124,220,82]/255,'MarkerSize',2);

            hold on;
            if(~bSkeleton)
                plot3(Up(:,2),Up(:,1),Up(:,3),'color',Colour);
                plot3(Down(:,2),Down(:,1),Down(:,3),'color',Colour);
                plot3(Left(:,2),Left(:,1),Left(:,3),'color',Colour);
                plot3(Right(:,2),Right(:,1),Right(:,3),'color',Colour);

                if(isfield(Branch,'UpLeft'))
                    plot3(UpLeft(:,2),UpLeft(:,1),UpLeft(:,3),'color',Colour);
                    plot3(DownRight(:,2),DownRight(:,1),DownRight(:,3),'color',Colour);
                    plot3(DownLeft(:,2),DownLeft(:,1),DownLeft(:,3),'color',Colour);
                    plot3(UpRight(:,2),UpRight(:,1),UpRight(:,3),'color',Colour);
                end

                for iP = 1:size(Boundary,1)

                    plot3(squeeze(Boundary(iP,2,:)),squeeze(Boundary(iP,1,:)),squeeze(Boundary(iP,3,:)),':','color',Colour);

                end
            end

        end
    end
end
axis equal;
hold off

end