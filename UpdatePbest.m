function Pcenter = UpdatePbest(Population)
% Update the ����λ��
   %% Select the Pcenter
    FrontNo=NDSort(Population.objs,1);
    temp=find(FrontNo==1);    %�ҵ��ڵ�һ�����и����λ��  
    temp1=Population.decs;
    temp2=temp1(temp,:); 
    Pcenter=(sum(temp2))/numel(temp);      %���������λ��(ʵ������Ӧֵ)

%     replace        = ~all(Pcenter2.objs>=Pcenter.objs,2);
%     Pcenter(replace) = Population(replace);
    
    
end