    function Population2=reProduce(Global,Population,Population1)
    
        position=Population.decs;
        fitness=Population.objs;
        
        position1=Population1.decs;
        
       %% parameter setting
       %  0.5(0)=0.7(0)
        pa=0.25; 
        %依概率Pa随机生成新的位置
        temp_flag3=rand(Global.N,Global.D)>pa;
%        temp_flag3=rand(Global.N,1)>pa;
        
%        new_position=position+rand*(position(randperm(Global.N),:)-position(randperm(Global.N),:));    


        FrontNo    = NDSort(Population.objs,inf);
        CrowdDis   = CrowdingDistance(Population.objs,FrontNo);
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        ePOP=Population.decs(MatingPool);
        temp1=rand(Global.N, 1);
        temp=repmat(temp1,1,Global.D);
%         new_position=position+rand*(ePOP(randperm(Global.N),:)-ePOP(randperm(Global.N),:)); 


%         temp2=rand(Global.N, Global.D)>0.5/Global.M;
%         temp2=rand(Global.N, Global.D)>1/Global.M;
         
        temp2=rand(Global.N, Global.D)>1/Global.M;
         
         
%         temp2=rand(Global.N, Global.D);

        % 2D   0.4(4)=0.3(4)=0.1(4)=0.7(4)=0.8(4)=0.2(5)<0.6(5)=0.5(5)
        % 4D   0.5(1)=0.1(1)=0.2(1)=0.6(1)=0.8(1)=0.4(1)=<0.25(2)=0.3(2) 
        % 6D   0.167(0)        
        
%         new_position=randi([0,1],1,Global.D).*(temp.*ePOP(randperm(Global.N),:)+(1-temp).*ePOP(randperm(Global.N),:));
%         new_position=temp2.*(temp.*ePOP(randperm(Global.N),:)+(1-temp).*ePOP(randperm(Global.N),:));


          new_position=ePOP+temp2.*temp.*(position1(randperm(Global.N),:)-position1(randperm(Global.N),:));
          
          
          
 
%        stepsize=0.01*step.*(randi([-1,1],1,Global.D).*(Global.upper-Global.lower)*(1-Global.evaluated/Global.evaluation));        
%         new_position=modify_position(new_position,xmax,xmin);
%         new_fitness=calculation_fitness(new_position);
       [position2,fitness1,~] = Global.problem('value',Global,new_position);
     
%         temp_flag4=new_fitness<fitness;
        temp_flag4=Domination(Global.N,Global.M,fitness1,fitness);
        temp_flag5=temp_flag3&repmat(temp_flag4,1,Global.D);
%         temp_flag6=repmat(temp_flag5,1,Global.D); 
        position=position2.*temp_flag5+position.*(~temp_flag5);
%         position=position2.*temp_flag6+position.*(~temp_flag6);
%         fitness=new_fitness.*temp_flag5+fitness.*(~temp_flag5);
        Population2=INDIVIDUAL(position);
    end