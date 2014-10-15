function [Vessels] = VesselAnalysis3D(Skeleton, VesselSeg)
%   Function to open a segmentation or skeleton of vasculature and
%   decompose this into its network structure. Outputted is a hierarchical
%   tree along with various descriptors of the tree. 


%   Varies whether you will read in a vessel segmentation or a vessel
%   skeleton. 
bSkeleton =0;
OutputName = 'Vessels';
folder1 = cd;

FrangiThreshold = 0.2;

if (nargin<1)
    if (bSkeleton == 1)
        [fname, folder1] = uigetfile('*.nii', 'Open Skeleton Image');
        Skeleton = load_untouch_nii([folder1, fname]);
        [fname, folder1] = uigetfile([folder1,'*.nii'], 'Open Vessel Segmentation');
        VesselSeg = load_untouch_nii([folder1, fname]);
        VesselSegImg = padarray(VesselSeg.img,[1,1,1]);
        SkeletonImg = Skeleton.img;
        SkeletonImg = padarray(SkeletonImg,[1,1,1]);
    else
        [fname, folder1] = uigetfile('*.nii', 'Open Vessel Segmentation');
        VesselSeg = load_untouch_nii([folder1, fname]);

        SkeletonImg = Skeleton3D(VesselSeg.img > 0.5);
        Skeleton = VesselSeg;
        Skeleton.img = SkeletonImg;
        SkeletonImg = padarray(SkeletonImg,[1,1,1]);
        VesselSegImg = padarray(VesselSeg.img,[1,1,1]);
        save_untouch_nii(Skeleton,[folder1, 'SkeletonImage.nii']);
    end
end

if(nargin == 1)
  
    VesselSegImg    = double(Skeleton.img);
    mins            = min(VesselSegImg(:));
    maxs            = max(VesselSegImg(:));
    
    SkeletonImg     = Skeleton3D((VesselSegImg-mins)./(maxs-mins) > FrangiThreshold);
    SkeletonImg     = padarray(SkeletonImg,[1,1,1]);
    VesselSegImg    = padarray(VesselSegImg,[1,1,1]);
end

if (nargin == 2)
   
    VesselSegImg = padarray(VesselSeg,[1,1,1]);
    SkeletonImg = padarray(Skeleton,[1,1,1]);
    
    
end

%Connected components are labeled. 
[SkelLabel nLabels] = bwlabeln(SkeletonImg,26);
if(isfield(Skeleton, 'hdr'))
    PixDims = Skeleton.hdr.dime.pixdim(2:4);
else
    PixDims = [1,1,1];
end

clear VesselSeg Skeleton;

ViewSkeleton(SkeletonImg)

TrueLabels = [];
nTrueLabels = 0;
Volume = 0;

for iL = 1:nLabels  
    
    if (nnz(SkelLabel == iL) > 50)
        TrueLabels = [TrueLabels iL];
        nTrueLabels = nTrueLabels + 1;
        SingleVessel = SkelLabel == iL;
        
        nPixels = nnz(SingleVessel);
        Volume = Volume + nPixels*prod(PixDims);
        
        %   Information on the connectedness of each network. Branchpoints,
        %   Endpoints etc. 
        [VesselData] = AnalyseVessel(SingleVessel);
        Vessels{nTrueLabels} = VesselData;
        Vessels{nTrueLabels}.nPixels = nPixels;
        Vessels{nTrueLabels}.Volume = Volume;
    end
    
end

for iL = 1:nTrueLabels
   
    SingleVessel = SkelLabel == TrueLabels(iL);
    %   Full tree structure of vessels is described along with descriptors of the
    %   network
    [BranchData] = TrackBranching(SingleVessel,Vessels{iL});
    Vessels{iL}.Branching  =  BranchData;
    
end

%   Local thickness of the vessels at each skeleton point can be calculated
%   by considering vectors fired out at various angles locally normal to the
%   skeleton. By considering the point of intersection between the vector
%   and a mesh element of the isosurface of the segmentation we can
%   calculate the thickness of the vessel in each of these directions. 
[Vessels] = CalculateThickness(Vessels,VesselSegImg)

PlotBranch(Vessels,1,1,1);
if(nargout == 0)
    save([folder1 OutputName '.mat'],'Vessels');
end

end



function [VesselData] = AnalyseVessel(SingleVessel)
%   Function takes the skeleton of a single vessel within the segmentation
%   and identifies the degree of connectivity for each point on the branch.
%   1 neighbour = end point
%   2 neighbours = regular point
%   3 neightbours = branching point
%   Returns number and location of branching/end points


ImSize = size(SingleVessel);

Directions = npermutek([-1,0,1],3);
DirectionSize = size(Directions,1);

[Points] = find(SingleVessel == 1);
Points = unique(Points);

[x,y,z] = ind2sub(ImSize,Points);
Subs = [x,y,z];

for iP = 1:size(Subs,1)
   
    Point = Subs(iP,:);
    
    Neighbours = repmat(Point,[DirectionSize,1]) + Directions;
    NeighbourInd = sub2ind(ImSize,Neighbours(:,1),Neighbours(:,2),Neighbours(:,3));
    
    NumNeighbours(iP) = nnz(SingleVessel(NeighbourInd)) - 1;
    
end

EndPoints = find(NumNeighbours == 1);

VesselData.EndPoints.Idx = EndPoints;
VesselData.EndPoints.Subs = Subs(EndPoints,:);

BranchPoints = find(NumNeighbours > 2);

VesselData.BranchPoints.Idx = BranchPoints;
VesselData.BranchPoints.Subs = Subs(BranchPoints,:);

end

function [BranchingData] = TrackBranching(SingleVessel,Vessels)
%   Function traces the skeleton of a single vessel along from its base
%   points explicitly collecting the points making up each branch and
%   saving the hierarchical structure of the tree. 

TmpVessel = SingleVessel;
ImSize = size(SingleVessel);
EndPoints = Vessels.EndPoints.Subs;
BranchPoints = Vessels.BranchPoints.Subs;
nBranches = length(EndPoints) + length(BranchPoints);

Directions = npermutek([-1,0,1],3);
nDirections = size(Directions,1);

%--------------------------------------------------------------------------
%   Choose Base point - Currently considered to be the end point with the
%   smallest x-coordinate.
%--------------------------------------------------------------------------
[xMin xMinIdx] = min(EndPoints(:,1));
Base = EndPoints(xMinIdx,:);
%--------------------------------------------------------------------------

SingleVessel(sub2ind(ImSize,Base(1),Base(2),Base(3))) = 10;
BInds = sub2ind(ImSize,BranchPoints(:,1),BranchPoints(:,2),BranchPoints(:,3));
EInds = sub2ind(ImSize,EndPoints(:,1),EndPoints(:,2),EndPoints(:,3));
SingleVessel(BInds) = 2;
SingleVessel(EInds) = 3;


%   Algorithm iterates through the tree structure following a breadth-first
%   search. For each branch the points which make it up, its parent and its
%   children are saved. 

CurrLoc = Base;
BranchingData.nBranches = 1;
BranchingData.CurrentBranch = 1;
BranchingData.Branches{1}.Points = Base;
[BranchingData,SingleVessel] = TrackBranch(CurrLoc,SingleVessel,BranchingData);
BranchingData.CurrentBranch = BranchingData.CurrentBranch + 1;


while(BranchingData.CurrentBranch <= BranchingData.nBranches)
    
    CurrLoc = BranchingData.Branches{BranchingData.CurrentBranch}.Points;
    [BranchingData,SingleVessel] = TrackBranch(CurrLoc,SingleVessel,BranchingData);
    BranchingData.CurrentBranch = BranchingData.CurrentBranch + 1;
    
end




nPruned = 1;

while(nPruned>0)
%   Tree Pruning - small branches can remain as artifacts of the
%   skeletonization procedure rather than part of the true structure. These
%   are removed in pruning. This potentially leaves branches with a single
%   child, in this case the parent branch inherits the points and children
%   which the single child comprised of. 
[BranchingData nPruned] = PruneBranches(BranchingData);
end

for iB = 1:BranchingData.nBranches
    %   More indepth measures calculated for each branch using the
    %   skeleton e.g. tortuosity
    [Tortuosity,VesselLength,SmoothedPoints,ChordLengthRatio] = CalculateTortuosity(BranchingData.Branches{iB}.Points);
    BranchingData.Branches{iB}.Tortuosity = Tortuosity;
    BranchingData.Branches{iB}.VesselLength = VesselLength;
    BranchingData.Branches{iB}.SmoothedPoints = SmoothedPoints;
    BranchingData.Branches{iB}.ChordLengthRatio = ChordLengthRatio;
end

end

function[BranchingData,SingleVessel] = TrackBranch(Start, SingleVessel,BranchingData)

Directions = npermutek([-1,0,1],3);
Directions(ismember(Directions,[0,0,0],'rows'),:) = [];
nDirections = size(Directions,1);
ImSize = size(SingleVessel);
CurrentBranch = BranchingData.CurrentBranch;
Points = BranchingData.Branches{CurrentBranch}.Points;

Point = Start;

Neighbours = repmat(Point,[nDirections,1]) + Directions;
NeighbourInd = sub2ind(ImSize,Neighbours(:,1),Neighbours(:,2),Neighbours(:,3));
NeighbourVal = SingleVessel(NeighbourInd);

while(nnz(NeighbourVal) == 1)
    
    SingleVessel(Point(1),Point(2),Point(3)) = 0;
    
    Point = Neighbours(NeighbourVal,:);
    Points = [Points;Point];
    
    Neighbours = repmat(Point,[nDirections,1]) + Directions;
    NeighbourInd = sub2ind(ImSize,Neighbours(:,1),Neighbours(:,2),Neighbours(:,3));
    NeighbourVal = SingleVessel(NeighbourInd);
    
end

if(nnz(NeighbourVal) > 1)
    SingleVessel(Point(1),Point(2),Point(3)) = 0;
    nBranchesFound = nnz(NeighbourVal);
    BranchingData.Branches{CurrentBranch}.Points = [Points;Point];
    BranchBases = Neighbours(NeighbourVal,:);
    BranchingData.Branches{CurrentBranch}.Children = (1:nBranchesFound) + BranchingData.nBranches;
        for iB = 1:nBranchesFound
            NewBranchNum = iB + BranchingData.nBranches;
            BranchingData.Branches{NewBranchNum}.Points = BranchBases(iB,:);
            BranchingData.Branches{NewBranchNum}.Parent = CurrentBranch;
        end
    BranchingData.nBranches = BranchingData.nBranches + nBranchesFound;
end

if(nnz(NeighbourVal) == 0)
    BranchingData.Branches{CurrentBranch}.Points = [Points;Point];
    nBranchesFound = 0;
    BranchBases = [];
end

end

function [BranchingData,nPruned] = PruneBranches(BranchingData)

    PruneList = [];
    KeepList = [];
    
    nPruned = 0;
    
    for iB = 1:numel(BranchingData.Branches)
        
        if (size(BranchingData.Branches{iB}.Points,1) < 5 && ~isfield(BranchingData.Branches{iB},'Children'))
            PruneList = [PruneList, iB];
        else
            KeepList = [KeepList,iB];
        end
        
    end
    
    BranchingData.Branches(PruneList) = [];
    SwitchTable = zeros(1,max([KeepList,PruneList]));
    SwitchTable(KeepList) = 1:length(KeepList);
    nPruned = nPruned + length(PruneList);
    
    for iB = 1:numel(BranchingData.Branches)
       
        if (isfield(BranchingData.Branches{iB},'Parent'))
            Parent = BranchingData.Branches{iB}.Parent;
            BranchingData.Branches{iB}.Parent = SwitchTable(Parent);
        end
        
        if (isfield(BranchingData.Branches{iB},'Children'))
            Children = BranchingData.Branches{iB}.Children;
            NewChildren = SwitchTable(Children);
            NewChildren = NewChildren(NewChildren~=0);
            if numel(NewChildren) > 0
                BranchingData.Branches{iB}.Children = NewChildren;
            else
                BranchingData.Branches{iB} = rmfield(BranchingData.Branches{iB},'Children');
            end
        end
        
    end
    

    BranchingData.nBranches = length(KeepList);
    
    if isfield(BranchingData,'CurrentBranch')
        BranchingData = rmfield(BranchingData,'CurrentBranch');
    end
    
    PruneList = [];
    KeepList = [];
    
    for iB = BranchingData.nBranches:-1:1
        
        BranchData = BranchingData.Branches{iB};
        
        if (isfield(BranchData,'Children')&&numel(BranchData.Children) == 1)
            
            ChildIdx = BranchData.Children;
            
            PruneList = [PruneList ChildIdx];
            
            BranchData.Points = [BranchData.Points; BranchingData.Branches{ChildIdx}.Points];
            
            if isfield(BranchingData.Branches{ChildIdx},'Children')
                Children = BranchingData.Branches{ChildIdx}.Children;
                BranchData.Children = Children;
                for iC = 1:length(Children)
                    BranchingData.Branches{Children(iC)}.Parent = iB; 
                end
                    
            else
                BranchData = rmfield(BranchData,'Children');
            end
       
            
        end
        
        BranchingData.Branches{iB} = BranchData;
            
    end
    
    KeepList = setdiff(1:BranchingData.nBranches,PruneList);
    
    BranchingData.Branches(PruneList) = [];
    SwitchTable = zeros(1,max([KeepList,PruneList]));
    SwitchTable(KeepList) = 1:length(KeepList);
    BranchingData.Connectivity = zeros(length(KeepList));
    nPruned = nPruned + length(PruneList);
    
    for iB = 1:length(KeepList)
       
        if (isfield(BranchingData.Branches{iB},'Parent'))
            Parent = BranchingData.Branches{iB}.Parent;
            BranchingData.Branches{iB}.Parent = SwitchTable(Parent);
        end
        
        if (isfield(BranchingData.Branches{iB},'Children'))
            Children = BranchingData.Branches{iB}.Children;
            NewChildren = SwitchTable(Children);
            NewChildren = NewChildren(NewChildren~=0);
            if numel(NewChildren) > 0
                BranchingData.Branches{iB}.Children = NewChildren;
            else
                BranchingData.Branches{iB} = rmfield(BranchingData.Branches{iB},'Children');
            end
        end
        
        if isfield(BranchingData.Branches{iB},'Children')
            BranchingData.Connectivity(iB,BranchingData.Branches{iB}.Children) = 1;
        end
        
    end
    
    BranchingData.nBranches = length(KeepList);
        

end
    




