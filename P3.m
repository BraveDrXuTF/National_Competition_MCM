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
%%%%%%%%%%%%%%%%%%%%%%按适应度升序排列%%%%%%%%%%%%%%%%%%%%%%%
for np=1:NP
    [S_left(np),T_3{np},v{np},b{np},t_region{np},T_region{np},k_interval{np}]=S_compute(f(:,np),betaall);
end
[SortS_left,Index]=sort(S_left);                            
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
          [NS_left(np),NT_3{np},Nv{np},Nb{np},Nt_region{np},NT_region{np},Nk_interval{np}]=S_compute(nf(:,np),betaall);   
    end
    [NSortS_left,Index]=sort(NS_left);           
    NSortf=nf(:,Index);
    %%%%%%%%%%%%%%%%%%%%%%%%%产生新种群%%%%%%%%%%%%%%%%%%%%%%%%%%
    f1=[Sortf,NSortf];                %子代和父代合并
    S_left1=[SortS_left,NSortS_left];       %子代和父代的适应度值合并
    [SortS_left1,Index]=sort(S_left1);    %适应度按升序排列
    Sortf1=f1(:,Index);               %按适应度排列个体
    SortS_left=SortS_left1(1:NP);         %取前NP个适应度值
    Sortf=Sortf1(:,1:NP);             %取前NP个个体
    trace(gen)=SortS_left(1);           %历代最优适应度值
end
Bestf=Sortf(:,1)                     %最优个体 
trace(end)                            %最优值
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

%%%%%%%%%%%%%%%%%%画出最优的曲线%%%%%%%%%%%%

[S,T_3,v,b,t_region,T_region,k_interval]=S_compute(Bestf,betaall);
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

plot(P1_x,P1_y);






