function [ Tortuosity, VesselLength, SmoothedPoints, ChordLengthRatio ] = CalculateTortuosity( VesselPoints )
%Calculates the tortuosity of a given blood vessel by calculating the
%integral of the squared curvature parameterised by the curve length,
%divided by the total curve length. 

nPoints = size(VesselPoints,1);

if (nPoints > 2);
    SmoothedPoints = SmoothPoints(VesselPoints);
    VesselPoints = SmoothPoints(SmoothedPoints);
else
    SmoothedPoints = VesselPoints;
end

if (nPoints > 5)
    Dt = VesselPoints(2:end,:) - VesselPoints(1:end-1,:);
    DtNorm = sqrt(diag(Dt*Dt'));
    VesselLength = sum(DtNorm);
    RA = VesselPoints(1:end-2,:);
    RO = VesselPoints(2:end-1,:);
    RB = VesselPoints(3:end,:);
    
    A = repmat(DtNorm(1:end-1),[1,3]);
    B = repmat(DtNorm(2:end),[1,3]);

    
    DpDt =  (RB - RO)./(B)+...
            (RO - RA)./(A)-...
            (RB - RA)./(A + B);

    D2pDt2 =    2*RA./(A.*(A + B)) -...
                2*RO./(A.*B)+...
                2*RB./((A+B).*B);

    N1 = D2pDt2(:,3).*DpDt(:,2) - D2pDt2(:,2).*DpDt(:,3);
    N2 = D2pDt2(:,1).*DpDt(:,3) - D2pDt2(:,3).*DpDt(:,1);
    N3 = D2pDt2(:,2).*DpDt(:,1) - D2pDt2(:,1).*DpDt(:,2);

    D1 = DpDt(:,1).^2 + DpDt(:,2).^2 + DpDt(:,3).^2;

    K = sqrt(N1.^2 + N2.^2 + N3.^2)./D1.^(3/2);

    CurveSpacing = cumsum(DtNorm(2:end));
    Tortuosity = trapz(CurveSpacing(2:end-1),K(2:end-1).^2)./VesselLength;
    ChordLengthRatio = norm(SmoothedPoints(1,:) - SmoothedPoints(end,:))/VesselLength;
else 
    Tortuosity  = 0;
    VesselLength = nPoints;
    ChordLengthRatio = 1;
end

end
