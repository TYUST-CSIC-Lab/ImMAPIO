    function Population2=reProduce2(Global,Population,Population1)
    % �ֲ�����
        position=Population.decs;
        speed=Population.adds;
        fitness=Population.objs;       
        position1=Population1.decs;
        
       %%�����һ֧��������λ��
        FrontNo    = NDSort(Population.objs,1);
        CrowdDis   = CrowdingDistance(Population.objs,FrontNo);
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        ePOP=Population.decs(MatingPool);
        PositionCenter=(sum(sum(ePOP)))/(numel(ePOP));
        new_position=PositionCenter-position1(randperm(Global.N),:);       
% �Ը����ÿ��ά�Ƚ��и���
       [position2,fitness1,~] = Global.problem('value',Global,new_position);
        temp_flag4=Domination(Global.N,Global.M,fitness1,fitness);
        temp_flag5=repmat(temp_flag4,1,Global.D);
        position=position2.*temp_flag5+position.*(~temp_flag5);
        Population2=INDIVIDUAL(position,speed);
    end