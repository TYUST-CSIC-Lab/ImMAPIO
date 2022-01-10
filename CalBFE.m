function BFE = CalBFE(PopObj,Zmin)
% Calculate the FC value of each solution

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".
%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group

    [N,M]  = size(PopObj);
    fmin = repmat(min(PopObj,[],1),N,1);
    fmax = repmat(max(PopObj,[],1),N,1);
    PopObj = (PopObj - fmin)./(fmax-fmin);
    
    
        %% Calculate Cv
    dis = sqrt(sum(PopObj.^2,2));
    % Use 1-dis instead of 1-dis/sqrt(M)
    Cv  = 1 - dis;
    
    %% Calculate d1 and d2
    Cosine = 1 - pdist2(PopObj,ones(1,M),'cosine');
    d1     = dis.*Cosine;
    d2     = dis.*sqrt(1-Cosine.^2);
        %% Calculate the shifted distance between each two solutions
    sde = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            sde(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    
    SDE = min(sde,[],2);
    Cd  = (SDE-min(SDE))./(max(SDE)-min(SDE));
    
    %% Determine the value of alpha and beta of each solution
    alpha   = zeros(length(Cv),1);
    beta    = zeros(length(Cv),1);
    meanCd  = mean(Cd);
    meanCv  = mean(Cv);
    meand1  = mean(d1);
    meand2  = mean(d2);
    case111 = Cv >  meanCv & d1 <= meand1 & Cd <= meanCd;
    case112 = Cv >  meanCv & d1 <= meand1 & Cd >  meanCd;
    case121 = Cv >  meanCv & d1 >  meand1 & Cd <= meanCd;
    case122 = Cv >  meanCv & d1 >  meand1 & Cd >  meanCd;
    case211 = Cv <= meanCv & d1 <= meand1 & d2 >  meand2 & Cd <= meanCd;
    case212 = Cv <= meanCv & d1 <= meand1 & d2 >  meand2 & Cd >  meanCd;
    case221 = Cv <= meanCv &(d1 >  meand1 | d2 <= meand2)& Cd <= meanCd;
    case222 = Cv <= meanCv &(d1 >  meand1 | d2 <= meand2)& Cd >  meanCd;
    alpha(case111) = rand(sum(case111),1)*0.3+0.8; beta(case111) = 1;
    alpha(case112) = 1;   beta(case112) = 1;
    alpha(case121) = 0.6; beta(case121) = 1;
    alpha(case122) = 0.9; beta(case122) = 1;
    alpha(case211) = rand(sum(case211),1)*0.3+0.8; beta(case211) = rand(sum(case211),1)*0.3+0.8;
    alpha(case212) = 1;   beta(case212) = 1;
    alpha(case221) = 0.2; beta(case221) = 0.2;
    alpha(case222) = 1;   beta(case222) = 0.2;

    %% The BFE value of each solution
    BFE = alpha.*Cd + beta.*Cv;
    
%     %% Favorable weight
%     w     = zeros(N,M);
%     bound = any(PopObj==repmat(Zmin,N,1),2);
%     w(repmat(bound,1,M) & PopObj==0) = 1;
%     w(repmat(bound,1,M) & PopObj~=0) = 0;
%     w(~bound,:) = 1./PopObj(~bound,:)./repmat(sum(1./PopObj(~bound,:),2),1,M);
%     
%     %% Calculate the FC value
%     FC = max(w.*PopObj,[],2);
%     FC = max(FC,1e-6);
end