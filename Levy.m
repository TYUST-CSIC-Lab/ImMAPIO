    function Population2=Levy(Global,Population)
        tr=1;          %tr为过渡因子
        R=0.3;             %R为地图因子
%         Population= Population(min(1:ceil(end/3)*3,end));
        Population= Population([1:end,1:ceil(end/3)*3-end]);
%         fitness=Population.objs;=
        position1=Population.decs;       %位置
        [N,D]     = size(position1);
        speed1 = Population.adds(zeros(N,D));   %速度
%         c1=2.05;                    % Acceleration constants, usually between 0~2, c1 towards the best individual postion
%         c2=2.05;                       % Acceleration constants towards the best global position
%         phi=c1+c2;                 % Constriction
%         kapa=2/abs(2-phi-sqrt(phi^2-4*phi));      %   Constriction factor
        position2   = position1(1:N/3,:);
        speed2 = speed1(1:N/3,:);
        PositionCenter= position1(N/3+1:N/3*2,:);
        PositionBest  = position1(N/3*2+1:end,:);
%         %开始执行Levy飞行
%         %% 全局搜索更新     
%         beta=3/2;
%         sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
%         for k=1:Global.N
%             
%         %这是实现Levy飞行的一个简单方式，即对于随机步长来说，采用步长为1
%          u=randn(size(position1(k,:)))*sigma;
%          v=randn(size(position1(k,:)));
%          step=u./abs(v).^(1/beta);               %计算的step为服从Levy飞行的随机数
% %            stepsize=0.01*step.*(randi([-1,1],1,Global.D).*(Global.upper-Global.lower)*(1-(log(Global.evaluated)/log(Global.evaluation))));
%            stepsize=tr*step.*(randi([-1,1],1,Global.D).*(Global.upper-Global.lower)*(1-Global.evaluated/Global.evaluation));   %(randi([-1,1],1,Global.D)表示服从高斯分布的随机数
%          
%        %% 局部搜索更新    
%        %计算第一支配层的中心位置
%         FrontNo= NDSort(fitness,1);  %发现第一层前沿面，和最大前沿面数
%         CrowdDis   = CrowdingDistance(fitness,FrontNo);
%         MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
%         ePOP=position1(MatingPool);    %ePOP表示第一支配层所有个体的位置
%         PositionCenter=(sum(sum(ePOP)))/(numel(ePOP));      
%         %% 全局和局部实际的速度和位置更新
%          speed2(k,:)= speed1(k,:)+stepsize.*randn(size(position1(k,:)))+(PositionCenter-position1(k,:)).*(randi([-1,1],1,Global.D))*(Global.evaluated/Global.evaluation)*tr;                  %此处速度未加局部更新
%          Bestfitness=FrontNo(randperm(length(FrontNo),1)); 
%          position2(k,:)=position1(k,:)*(1-exp(-R*(abs(Bestfitness-fitness(k)))))+kapa*speed2(k,:); 
% %          position2(k,:)=position1(k,:)*(1-exp(-R*Global.evaluated))+speed2(k,:);   
%     
%         end
%        %% 实现各个维度的更新
%         [Decs2,Objs2,~] = Global.problem('value',Global,position2);
%         [Decs1,Objs1,~] = Global.problem('value',Global,position1);
%               
%         Flag=Domination(Global.N,Global.M,Objs2,Objs1);
%         new_Decs=Decs2.*repmat(Flag,1,Global.D)+Decs1.*repmat(~Flag,1,Global.D);
%         new_Speed= speed2;
%         Population2=INDIVIDUAL(new_Decs,new_Speed);  
    
        r1 = repmat(rand(N/3,1),1,D);
        r2 = repmat(rand(N/3,1),1,D);
        r3 = repmat(rand(N/3,1),1,D);
        C1 = repmat(rand(N/3,1)+1.5,1,D)>1/Global.M;
        C2 = repmat(rand(N/3,1)+1.5,1,D)>1/Global.M;
        C3 = repmat(rand(N/3,1)+1.5,1,D)>1/Global.M;
% (1-(log(Global.evaluated)/log(Global.evaluation))
        speed3= speed2*(exp(-R*Global.evaluated))+tr*(1-(log(Global.evaluated)/log(Global.evaluation))).*r1.*(PositionBest-position2).*C1+tr*(log(Global.evaluated)/log(Global.evaluation)).*r2.*(PositionCenter-position2).*C2+C3.*r3.*(PositionBest-PositionCenter);                  %此处速度未加局部更新
        position3=position2+speed3;
%         Population2= INDIVIDUAL(position3,speed3);  
       %% 实现各个维度的更新
        [Decs2,Objs2,~] = Global.problem('value',Global,position3);
        [Decs1,Objs1,~] = Global.problem('value',Global,position2);

        Flag=Domination(N/3,Global.M,Objs2,Objs1);
        new_Decs=Decs2.*repmat(Flag,1,Global.D)+Decs1.*repmat(~Flag,1,Global.D);
        new_Speed= speed3;
        Population2=INDIVIDUAL(new_Decs,new_Speed); 
        
    end 