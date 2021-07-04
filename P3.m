clc,clear
load betaall.mat
D=5;                                            %������Ŀ��Ԫ��
NP=100;                                         %Ⱦɫ����Ŀ
G=1000;                                         %����Ŵ�����
f=zeros(D,NP);                                  %��ʼ��Ⱥ���ռ�
nf=zeros(D,NP);                                  %����Ⱥ���ռ�
Pc=0.8;                                         %�������
Pm=0.1;                                         %�������
f=(rand(D,NP)-0.5)*2;                         %�����ó�ʼ��Ⱥ-1,1
%%%%%%%%%%%%%%%%%%%%%%����Ӧ����������%%%%%%%%%%%%%%%%%%%%%%%
for np=1:NP
    [S_left(np),T_3{np},v{np},b{np},t_region{np},T_region{np},k_interval{np}]=S_compute(f(:,np),betaall);
end
[SortS_left,Index]=sort(S_left);                            
Sortf=f(:,Index);
%%%%%%%%%%%%%%%%%%%%%%%�Ŵ��㷨ѭ��%%%%%%%%%%%%%%%%%%%%%%%%%%
for gen=1:G
    %%%%%%%%%%%%%%���þ�����������ѡ�񽻲����%%%%%%%%%%%%%%%%
    Emper=Sortf(:,1);                      %����Ⱦɫ��
    NoPoint=round(D*Pc);                   %ÿ�ν����ĸ���
    PoPoint=randi(D,NoPoint,NP/2);   %��������λ��
    nf=Sortf;
    for i=1:NP/2
        nf(:,2*i-1)=Emper;
        nf(:,2*i)=Sortf(:,2*i);
        for k=1:NoPoint
            nf(PoPoint(k,i),2*i-1)=nf(PoPoint(k,i),2*i);
            nf(PoPoint(k,i),2*i)=Emper(PoPoint(k,i));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%�������%%%%%%%%%%%%%%%%%%%%%%%%%
    for m=1:NP
        for n=1:D
            r=rand(1,1);
            if r<Pm
                nf(n,m)=2*(rand(1,1)-0.5);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%����Ⱥ����Ӧ����������%%%%%%%%%%%%%%%%%%
    for np=1:NP 
          [NS_left(np),NT_3{np},Nv{np},Nb{np},Nt_region{np},NT_region{np},Nk_interval{np}]=S_compute(nf(:,np),betaall);   
    end
    [NSortS_left,Index]=sort(NS_left);           
    NSortf=nf(:,Index);
    %%%%%%%%%%%%%%%%%%%%%%%%%��������Ⱥ%%%%%%%%%%%%%%%%%%%%%%%%%%
    f1=[Sortf,NSortf];                %�Ӵ��͸����ϲ�
    S_left1=[SortS_left,NSortS_left];       %�Ӵ��͸�������Ӧ��ֵ�ϲ�
    [SortS_left1,Index]=sort(S_left1);    %��Ӧ�Ȱ���������
    Sortf1=f1(:,Index);               %����Ӧ�����и���
    SortS_left=SortS_left1(1:NP);         %ȡǰNP����Ӧ��ֵ
    Sortf=Sortf1(:,1:NP);             %ȡǰNP������
    trace(gen)=SortS_left(1);           %����������Ӧ��ֵ
end
Bestf=Sortf(:,1)                     %���Ÿ��� 
trace(end)                            %����ֵ
figure
plot(trace)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')

%%%%%%%%%%%%%%%%%%�������ŵ�����%%%%%%%%%%%%

[S,T_3,v,b,t_region,T_region,k_interval]=S_compute(Bestf,betaall);
%%%%%%%%%%%��ͼ%%%%%%%%%%%%%%%%%%%%%
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






