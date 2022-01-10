function ImMAPIO(Global)
% <algorithm> <O-Z>
% A Reference Vector Guided Evolutionary Algorithm for Many-objective
% Optimization
% alpha ---   2 --- The parameter controlling the rate of change of penalty
% fr    --- 0.1 --- The frequency of employing reference vector adaptation

%--------------------------------------------------------------------------
% The copyright of the PlatEMO belongs to the BIMK Group. You are free to
% use the PlatEMO for research purposes. All publications which use this
% platform or any code in the platform should acknowledge the use of
% "PlatEMO" and reference "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu
% Jin, PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective
% Optimization, 2016".

%--------------------------------------------------------------------------

% Copyright (c) 2016-2017 BIMK Group
    
    %t = Global.ParameterSet(0);
   %% Parameter setting
   [alpha,delta] = Global.ParameterSet(0.4,0.1);
   nSample = Global.ParameterSet(10000);
  
   rate = Global.ParameterSet(0.5);
   %% Generate random population
   [Z,Global.N] = UniformPoint(Global.N,Global.M);
    Population  = Global.Initialization();   %% elitite pop   父代
%     Population1 = Global.Initialization();   %% current pop子代种群
    
    %HyPE
    RefPoint = zeros(1,Global.M) + max(Population.objs)*1.2;
     
    %KnEA
    FrontNo    = NDSort(Population.objs,inf);
    KneePoints = zeros(1,Global.N);     % Set of knee points
    r          = -ones(1,2*Global.N);	% Ratio of size of neighorhood
    t          = -ones(1,2*Global.N);	% Ratio of knee points
   
    
    Zmin1       = min(Population.objs,[],1);      
    Fitness    = CalFitness(Population.objs);
    R  = GenerateRefPoints(Population,delta*(max(Population.objs,[],1)-min(Population.objs,[],1)),alpha,Global.N);
    Archive    = UpdateArchive(Population(NDSort(Population.objs,1)==1),[],Global.N);
   
    %% Optimization优化初始种群
    %while Global.NotTermination(Population)  
       
        Population1=Population;
        Population4=Population;
        Population5=Population;
    %end 
    
    Pcenter      = UpdatePbest(Population);
    %FrontNo=NDSort(Population.objs,1);
    %temp=find(FrontNo==1);    %找到在第一层所有个体的位置  
    %temp1=Population.decs;
    %temp2=temp1(temp,:); 
    %Pcenter=(sum(temp2))/numel(temp);      %计算出中心位置(实际是适应值)

    Archive    = UpdateArchive(Population(NDSort(Population.objs,1)==1),[],Global.N);

    %% Optimization
    while Global.NotTermination(Archive)
%         Global.NC=ceil(Global.N-Global.N*(0.4^( Global.evaluated)));         %每次迭代后剩余的种群数目
        %选择
           
        %HyPE
        MatingPool1 = TournamentSelection(2,Global.N,-CalHV(Population1.objs,RefPoint,Global.N,nSample));
        Offspring1  = Global.Variation(Population1(MatingPool1));    
        Population1 = EnvironmentalSelection1([Population1,Offspring1],Global.N,RefPoint,nSample);
        Archive    = UpdateArchive(Archive,Population1,Global.N);
       
        %KnEA
        MatingPool4 = MatingSelection4(Population4.objs,FrontNo,KneePoints);
        Offspring4  = Global.Variation(Population4(MatingPool4));
        Population4 = [Population4,Offspring4];
        [FrontNo,MaxFNo]                = NDSort(Population4.objs,Global.N);
        [KneePoints,Distance,r,t]       = FindKneePoints(Population4.objs,FrontNo,MaxFNo,r,t,rate);
        [Population4,FrontNo,KneePoints] = EnvironmentalSelection4(Population4,FrontNo,MaxFNo,KneePoints,Distance,Global.N);   
        Archive    = UpdateArchive(Archive,Population4,Global.N);
        
        %VaEA
        MatingPool5 = randi(Global.N,1,Global.N);
        Offspring5  = Global.Variation(Population5(MatingPool5));    
        Population5 = EnvironmentalSelection5([Population5,Offspring5],Global.N);
        Archive    = UpdateArchive(Archive,Population5,Global.N);
        
       
 
      
        %生成父代种群
        Population = Global.Variation([Population,Pcenter,Archive(randi(ceil(length(Archive)/10),1,Global.N))],Global.N,@Levy);
%         Population = Global.Variation(Population,Global.N,@Levy);
        %Pcenter      = UpdatePbest(Population);
        
        Archive    = UpdateArchive(Archive,Population,Global.N);
%         S          = Global.Variation(Archive(1:length(Archive)),length(Archive),@EAreal);
        S          = Global.Variation(Archive([1:length(Archive),randi(ceil(length(Archive)/2),1,length(Archive))]),length(Archive),@EAreal);
        Archive    = UpdateArchive(Archive,S,Global.N);      
    end

end






