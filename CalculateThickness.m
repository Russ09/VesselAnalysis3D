function [Vessels] = CalculateThickness(Vessels,VesselSegImage)


disp('Calculating Vessel Thickness');
OutputMesh = isosurface(VesselSegImage,0.1);
OutputMesh = reducepatch(OutputMesh,0.5);
OutputMesh = smoothpatch(OutputMesh);

nDirections = 8;

for iV = 1:numel(Vessels)
    
   CurrVes = Vessels{iV};

   [Out1,Out2] = surfacelandmarkRuss(OutputMesh,CurrVes,1,7,3,nDirections,[]);
       
   for iB = 1:numel(CurrVes.Branching.Branches)
       Branch = CurrVes.Branching.Branches{iB};
       nP = 0;
       for iP = 2:(size(Branch.Points,1)-1)
           Point = Branch.SmoothedPoints(iP,:);
           Point = repmat(Point,[nDirections,1]);
           
           if (nnz(Out1.inter{iB}{iP-1}) == 3*nDirections)

               nP = nP + 1;
               Intercepts = Out1.inter{iB}{iP-1};
               Intercepts = Intercepts(:,[2,1,3]);
               Diff = Point - Intercepts;
               Norms = sqrt(diag(Diff*Diff'));
               Branch.Thickness(nP) = mean(Norms);
               Branch.Up(nP,:) = Intercepts(1,:);
               Branch.Down(nP,:) = Intercepts(2,:);
               Branch.Left(nP,:) = Intercepts(3,:);
               Branch.Right(nP,:) = Intercepts(4,:);
               
               if nDirections == 8
                   Branch.UpLeft(nP,:) = Intercepts(5,:);
                   Branch.DownRight(nP,:) = Intercepts(6,:);
                   Branch.UpRight(nP,:) = Intercepts(7,:);
                   Branch.DownLeft(nP,:) = Intercepts(8,:);
               end
               
           end
           
       end
       
%        if(isfield(Branch,'Up'))
%            Branch.Up = SmoothPoints(Branch.Up);
%            Branch.Down = SmoothPoints(Branch.Down);
%            Branch.Left = SmoothPoints(Branch.Left);
%            Branch.Right = SmoothPoints(Branch.Right);
%            
%            if nDirections == 8
%                Branch.UpLeft = SmoothPoints(Branch.UpLeft);
%                Branch.DownRight = SmoothPoints(Branch.DownRight);
%                Branch.DownLeft = SmoothPoints(Branch.DownLeft);
%                Branch.UpRight = SmoothPoints(Branch.UpRight);
%            end
%                
%        end
       CurrVes.Branching.Branches{iB} = Branch;
   end
   Vessels{iV} = CurrVes;
    
end



end
