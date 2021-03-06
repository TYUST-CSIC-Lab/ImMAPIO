function MatingPool = MatingSelection4(PopObj,FrontNo,KneePoints)
% The mating selection of KnEA

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group

    N = size(PopObj,1);
    
    %% Calculate the weighted distance of each solution
    Dis   = pdist2(PopObj,PopObj);
    Dis(logical(eye(length(Dis)))) = inf;
    Dis   = sort(Dis,2);
	Crowd = sum(Dis(1:3,:).*repmat((3:-1:1)',1,N));

    %% Binary tournament selection
    MatingPool = TournamentSelection(2,N,FrontNo,-KneePoints,-Crowd);
end