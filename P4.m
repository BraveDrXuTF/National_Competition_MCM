clc,clear
load betaall.mat
D=5;                                            %基因数目五元组
NP=100;                                         %染色体数目
G=1000;                                         %最大遗传代数
f=zeros(D,NP);                                  %初始种群赋空间
nf=zeros(D,NP);                                  %子种群赋空间
Pc=0.8;                                         %交叉概率
Pm=0.1;                                         %变异概率
f=(rand(D,NP)-0.5)*2;                         %随机获得初始种群-1,1
rou=0.7;
%%%%%%%%%%%%%%%%%%%%%%按适应度升序排列%%%%%%%%%%%%%%%%%%%%%%%
for np=1:NP
    [S(np),S_right(np),T_3{np},v{np},b{np},t_region{np},T_region{np},k_interval{np}]=S_compute_P4(f(:,np),betaall);
    
end
num=1;
for np=1:NP
    
    if S(np)<inf
        S_valid(num)=S(np);
        if S_right(np)<inf
            S_right_valid(num)=S_right(np);
        end
        S_delta_valid(np)=abs(S_right_valid(num)-S_valid(num));
        num=num+1;
    else
        fitness(np)=inf;
    end
    
    
end
maxS=max(S_valid);
minS=min(S_valid);
maxS_delta=max(S_delta_valid);
minS_delta=min(S_delta_valid);
for np=1:NP
    
    if S(np)<inf
        
        fitness(np)=rou*(S(np)-minS)/(maxS-minS+0.001)+(1-rou)*(S_delta_valid(np)-minS_delta)/(maxS_delta-minS_delta+0.001);
        
    end
    
end



[Sortfitness,Index]=sort(fitness);                            
Sortf=f(:,Index);
%%%%%%%%%%%%%%%%%%%%%%%遗传算法循环%%%%%%%%%%%%%%%%%%%%%%%%%%
for gen=1:G
    %%%%%%%%%%%%%%采用君主方案进行选择交叉操作%%%%%%%%%%%%%%%%
    Emper=Sortf(:,1);                      %君主染色体
    NoPoint=round(D*Pc);                   %每次交叉点的个数
    PoPoint=randi(D,NoPoint,NP/2);   %交叉基因的位置
    nf=Sortf;
    for i=1:NP/2
        nf(:,2*i-1)=Emper;
        nf(:,2*i)=Sortf(:,2*i);
        for k=1:NoPoint
            nf(PoPoint(k,i),2*i-1)=nf(PoPoint(k,i),2*i);
            nf(PoPoint(k,i),2*i)=Emper(PoPoint(k,i));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%变异操作%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:NP
        for n=1:D
            r=rand(1,1);
            if r<Pm
                nf(n,m)=2*(rand(1,1)-0.5);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%子种群按适应度升序排列%%%%%%%%%%%%%%%%%%
    for np=1:NP 
          [NS(np),NS_right(np),NT_3{np},Nv{np},Nb{np},Nt_region{np},NT_region{np},Nk_interval{np}]=S_compute_P4(nf(:,np),betaall);   
    end
    
    %%%%%%%%%%%%%%%%%%%%计算归一化适应度%%%%%%%%%%%%%%%%%%%%%%%%%
    num=1;
    for np=1:NP

        if NS(np)<inf
            NS_valid(num)=NS(np);
            if NS_right(np)<inf
                NS_right_valid(num)=NS_right(np);
            end
            NS_delta_valid(np)=abs(NS_right_valid(num)-NS_valid(num));
            num=num+1;
        else
            Nfitness(np)=inf;
        end


    end
    maxNS=max(NS_valid);
    minNS=min(NS_valid);
    maxNS_delta=max(NS_delta_valid);
    minNS_delta=min(NS_delta_valid);
    for np=1:NP

        if NS(np)<inf

            Nfitness(np)=rou*(NS(np)-minNS)/(maxNS-minNS+0.001)+(1-rou)*(NS_delta_valid(np)-minNS_delta)/(maxNS_delta-minNS_delta+0.001);

        end

    end
    
    
    
    
    
    [NSortfitness,Index]=sort(Nfitness);           
    NSortf=nf(:,Index);
    %%%%%%%%%%%%%%%%%%%%%%%%%产生新种群%%%%%%%%%%%%%%%%%%%%%%%%%%
    f1=[Sortf,NSortf];                %子代和父代合并
    fitness1=[Sortfitness,NSortfitness];       %子代和父代的适应度值合并
    [Sortfitness1,Index]=sort(fitness1);    %适应度按升序排列
    Sortf1=f1(:,Index);               %按适应度排列个体
    Sortfitness=Sortfitness1(1:NP);         %取前NP个适应度值
    Sortf=Sortf1(:,1:NP);             %取前NP个个体
    trace(gen)=Sortfitness(1);           %历代最优适应度值
end
Bestf=Sortf(:,1)                     %最优个体 
trace(end)                            %最优值
h=figure
plot(trace)
xlim([0,50])
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')
saveas(h,'4.jpg')
%%%%%%%%%%%%%%%%%%画出最优的曲线%%%%%%%%%%%%

[S,S_right,T_3,v,b,t_region,T_region,k_interval]=S_compute_P4(Bestf,betaall);
%%%%%%%%%%%画图%%%%%%%%%%%%%%%%%%%%%
num=1;
b(1)=-30+T_3(1);%173-25
k_qian=b(1)*betaall(1,2);
for t=0:0.5:18
    P1_x(num) = t;
    P1_y(num) = 25+(30-25)/(18)*t;
    num=num+1;
end
for t = t_region(1,1):0.5:t_region(5,2)
    
    for r =1:5
        if t>t_region(r,1)&t<t_region(r,2)
            if r==1
                
                P1_x(num) = t;
                P1_y(num) = T_3(r)-b(r).*exp(-betaall(r,2)*(t-t_region(r,1)));
                num=num+1;
            elseif r>1&&r<=5
                ychu(r) = y0(r-1)+k_interval(r-1)*(t_region(r,1)-t_region(r-1,2));
                b(r)=T_3(r)-ychu(r);
                P1_x(num) = t;
                P1_y(num) = T_3(r)-b(r).*exp(-betaall(r,2)*(t-t_region(r,1)));
                num=num+1;
           
            end
        elseif t>t_region(r,2)&t<t_region(r+1,1)
            k_interval(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)));
            P1_x(num) = t;
            y0(r) = T_3(r)-b(r).*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)));
            P1_y(num) = y0(r)+k_interval(r)*(t-t_region(r,2));
            num=num+1;
        end
    end
    
end
P1_x, P1_y
h=figure
plot(P1_x,P1_y);
saveas(h,'4df.jpg')





