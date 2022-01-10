function Pcenter = UpdatePbest(Population)
% Update the 中心位置
   %% Select the Pcenter
    FrontNo=NDSort(Population.objs,1);
    temp=find(FrontNo==1);    %找到在第一层所有个体的位置  
    temp1=Population.decs;
    temp2=temp1(temp,:); 
    Pcenter=(sum(temp2))/numel(temp);      %计算出中心位置(实际是适应值)

%     replace        = ~all(Pcenter2.objs>=Pcenter.objs,2);
%     Pcenter(replace) = Population(replace);
    
    
end