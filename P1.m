clc,clear



filename='data.xlsx';
t_xlrange='A2:A710';
T_xlrange='B2:B710';
xlname='Sheet1';
t=xlsread(filename,xlname,t_xlrange);
T=xlsread(filename,xlname,T_xlrange);



s_region=[25,197.5;202.5,233;238,268.5;273.5,339.5;344.5,410.5];

beta0all=[30,0.1;30,0.1;30,0.1;30,0.1;-30,0.01];
% beta0all=[30,0.1,200;30,0.1,200;30,0.1,200;30,0.1,200;-30,0.01,25];
% beta0all=[30,0.1,0.1,0.1;30,0.1,0.1,0.1;30,0.1,0.1,0.1;30,0.1,0.1,0.1;-30,0.01,0.1,0.1];
betaall=[];
T_air_param_region=[175,195,235,255,25];

v=70/60;
t_region=s_region/v;



region_tdata={};
region_Tdata={};
for region = 1:5 
    region_tdata{region}=t(t>=t_region(region,1)&t<=t_region(region,2))-t_region(region,1);
    region_Tdata{region}=T(t>=t_region(region,1)&t<=t_region(region,2));
    
end


for region =1:5
    x=region_tdata{region};
    y=region_Tdata{region};
    beta0=beta0all(region,:)';
    myfunc = inline([num2str(T_air_param_region(region)),'-beta(1).*exp(-beta(2).*x)'],'beta','x');
%     myfunc = inline(beta(3),'-beta(1).*exp(-beta(2).*x)'],'beta','x');
%     myfunc = inline(['beta(3).*x.^2+beta(4).*x','+',num2str(T_air_param_region(region)),'-beta(1).*exp(-beta(2).*x)'],'beta','x');
    betaall(region,:) = nlinfit(x,y,myfunc,beta0)';
    
end

h=figure(1)
for region = 1:5
    
    x=region_tdata{region};
    y=region_Tdata{region};
    beta=betaall(region,:)';
    y1=T_air_param_region(region)-beta(1).*exp(-beta(2).*x);
    plot(x,y,'-r',x,y1,'-b');
    legend('真实炉温曲线','拟合炉温曲线')
    saveas(h,[num2str(region) '.png'])
end


T_1=[173,198,230,257,25];
v_1=78/60;
t_region_1=s_region/v_1;
k_interval = [];
b=[];%代替beta(1)
y0=[];%用于直线计算的
ychu=[];%用于曲线计算的
%%%%%%%%%%%%%%小温区3中点111.25cm,温度173，6中点217.75cm,温度198,7中点253.25，温度230%%%%%%%%%%%%%%%%%%%%

s_small_region = 25 : 35.5 : 25+35.5*11;
temp=s_small_region(2:end)-5;
s_small_regions(1,:)=s_small_region(1:end-1);
s_small_regions(2,:)=temp;
s_small_regions=s_small_regions';
z_3=mean(s_small_regions(3,:));
z_6=mean(s_small_regions(6,:));
z_7=mean(s_small_regions(7,:));
e_8=s_small_regions(8,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num=1;
b(1)=-30+T_1(1);%173-25
k_qian=b(1)*betaall(1,2);
for t=0:0.5:19.0
    P1_x(num) = t;
    P1_y(num) = 25+(30-25)/(19.5)*t;
    num=num+1;
end
for t = 19.5:0.5:315.5
    
    for r =1:5
        if t>t_region_1(r,1)&t<t_region_1(r,2)
            if r==1
                
                P1_x(num) = t;
                P1_y(num) = T_1(r)-b(r).*exp(-betaall(r,2)*(t-t_region_1(r,1)));
                num=num+1;
            elseif r>1&&r<=5
                ychu(r) = y0(r-1)+k_interval(r-1)*(t_region_1(r,1)-t_region_1(r-1,2));
                b(r)=T_1(r)-ychu(r);
                P1_x(num) = t;
                P1_y(num) = T_1(r)-b(r).*exp(-betaall(r,2)*(t-t_region_1(r,1)));
                num=num+1;
           
            end
        elseif t>t_region_1(r,2)&t<t_region_1(r+1,1)
            k_interval(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region_1(r,2)-t_region_1(r,1)));
            P1_x(num) = t;
            y0(r) = T_1(r)-b(r).*exp(-betaall(r,2)*(t_region_1(r,2)-t_region_1(r,1)));
            P1_y(num) = y0(r)+k_interval(r)*(t-t_region_1(r,2));P1_x, P1_y
            num=num+1;
        end
    end
    
end
P1_x, P1_y

plot(P1_x,P1_y);



T_3z=T_1(1)-b(1).*exp(-betaall(1,2)*(z_3-s_region(1,1))/v_1);
T_6z=T_1(2)-b(2).*exp(-betaall(2,2)*(z_6-s_region(2,1))/v_1);
T_7z=T_1(3)-b(3).*exp(-betaall(3,2)*(z_7-s_region(3,1))/v_1);
T_8e=T_1(4)-b(4).*exp(-betaall(4,2)*(e_8-s_region(4,1))/v_1);
for region = 1:5
    x=region_tdata{region};
    y=region_Tdata{region};
    beta=betaall(region,:)';
    y1=T_air_param_region(region)-beta(1).*exp(-beta(2).*x);
    plot(x,y,'-r',x,y1,'-b');
end
