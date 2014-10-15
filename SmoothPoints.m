function [ SmoothedPoints ] = SmoothPoints( VesselPoints )


V1 = VesselPoints(:,1);
V2 = VesselPoints(:,2);
V3 = VesselPoints(:,3);

VesselPoints = [V1(~isnan(V1)),V2(~isnan(V2)),V3(~isnan(V3))];

nPoints = size(VesselPoints,1);
SmoothedPoints(1,:) = VesselPoints(1,:);
SmoothedPoints(nPoints,:) = VesselPoints(nPoints,:);

if nPoints > 5
    SmoothedPoints(2,:) = (VesselPoints(1,:) + VesselPoints(2,:) + VesselPoints(3,:))/3; 
    SmoothedPoints(3:end-2,:) = (VesselPoints(1:end-4,:) + VesselPoints(2:end-3,:) +...
                                VesselPoints(3:end-2,:) + VesselPoints(4:end-1,:) +...
                                VesselPoints(5:end,:))/5;
    SmoothedPoints(nPoints-1,:) = (VesselPoints(nPoints-2,:) + VesselPoints(nPoints-1,:) + VesselPoints(nPoints,:))/3; 
else
    SmoothedPoints(2:end-1,:) = (VesselPoints(1:end-2,:) + VesselPoints(2:end-1,:) + VesselPoints(3:end,:))/3;
end

end

