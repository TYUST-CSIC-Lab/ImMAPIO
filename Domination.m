    function D_index=Domination(N,M,Objs1,Objs2)
    % N��ʾ��Ⱥ��С��M��ʾĿ�����
    D_index=false(N,1);
    for i=1:N       
        Mark=0;
        for j=1:M
            Objs1(i,j)<Objs2(i,j);
            Mark=Mark+1;
        end       
        if  Mark==M
            D_index(i,1)=true;
        end     
    end
    end