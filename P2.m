clc,clear
load betaall.mat
T_2=[182,203,237,254,25];

s_region=[25,197.5;202.5,233;238,268.5;273.5,339.5;344.5,410.5];
t_region={};
T_region=[];
k_interval=[];
b(1)=-30+T_2(1);%173-25
flag=1;
v_ok = [];
num=1;
for v = 65/60:0.001:100/60
    t_region=s_region/v;
    for r = 1:5
        if r==1            
            T_region(r,1)=30;
            T_region(r,2)=b(1)*(1-exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1))))+T_region(r,1);
            k_interval(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)));
            T_region(r+1,1) = T_region(r,2) + k_interval(r)*(t_region(r+1,1)-t_region(r,2));
            b(r+1) = T_2(r+1)-T_region(r+1,1);
        elseif r<5
            T_region(r,2)=b(r)*(1-exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1))))+T_region(r,1);
            k_interval(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)));
            T_region(r+1,1) = T_region(r,2) + k_interval(r)*(t_region(r+1,1)-t_region(r,2));
            b(r+1) = T_2(r+1)-T_region(r+1,1);
        elseif r==5
            T_region(r,2)=b(r)*(1-exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1))))+T_region(r,1);
        end
    end

%%%%%%%%%%%%%%ÅÐ¶ÏÐ±ÂÊ%%%%%%%%%%%%%%%%%%
    for r=1:5
    	k_max(r)=abs(b(r)*betaall(r,2));
    	if k_max(r)>3
    		flag = 0;
    	end
    end
    
%%%%%%%%%%%%%ÅÐ¶Ï150-190Ê±¼ä%%%%%%%%%%%%
	for r=1:4
		if 150>T_region(r,1)&150<T_region(r,2)
			t_150_190(1)=-1/betaall(r,2)*log((T_2(r)-150)/b(r))+t_region(r,1);
		
		elseif 150>T_region(r,2)&150<T_region(r+1,1)			
			t_150_190(1)=(150-T_region(r,2))/k_interval(r)+t_region(r,2);
		end	
		
		if 190>T_region(r,1)&190<T_region(r,2)
			t_150_190(2)=-1/betaall(r,2)*log((T_2(r)-190)/b(r))+t_region(r,1);
		
		elseif 190>T_region(r,2)&190<T_region(r+1,1)			
			t_150_190(2)=(190-T_region(r,2))/k_interval(r)+t_region(r,2);
		end	
    end
    if ~(60<t_150_190(2)-t_150_190(1)&120>t_150_190(2)-t_150_190(1))
    	flag = 0;
    end
%%%%%%%%%%%%%ÅÐ¶Ï217Ê±¼ä%%%%%%%%%%%%%%
	for r=1:4
		if 217>T_region(r,1)&217<T_region(r,2)
			t_217(1)=-1/betaall(r,2)*log((T_2(r)-217)/b(r))+t_region(r,1);
		
		elseif 217>T_region(r,2)&217<T_region(r+1,1)			
			t_217(1)=(217-T_region(r,2))/k_interval(r)+t_region(r,2);
		end	
	end
	if 217>T_region(5,2)&217<T_region(5,1)
		t_217(2)=-1/betaall(5,2)*log((T_2(5)-217)/b(5))+t_region(5,1);
	end
    if ~(t_217(2)-t_217(1)>40&t_217(2)-t_217(1)<90)
    	flag = 0;
    end
%%%%%%%%%%%%ÅÐ¶Ï·åÖµ%%%%%%%%%%%%%%%%%%
	T_max=T_region(5,1);
	if ~(T_max>240&T_max<250)
		flag=0;
	end

	if flag==1
		v_ok(num)=v;
        b_ok{num}=b;
        t_region_ok{num}=t_region;
        k_interval_ok{num}=k_interval;
        t_217_ok{num} = t_217;
        t_150_190_ok{num} = t_150_190;
        T_max_ok{num}=T_max;
		num=num+1;
    end


    flag=1;

end
v_max=0;
[v_max,ind]=max(v_ok);
b_vm=b_ok{ind};
t_region_vm=t_region_ok{ind};
k_interval_vm=k_interval_ok{ind};
t_217_vm = t_217_ok{ind};
t_150_190_vm = t_150_190_ok{ind};
T_max_ok_vm = T_max_ok{ind};

%%%%%%%%%%%»­Í¼%%%%%%%%%%%%%%%%%%%%%
num=1;
b(1)=-30+T_2(1);%173-25
k_qian=b(1)*betaall(1,2);
for t=0:0.5:18
    P1_x(num) = t;
    P1_y(num) = 25+(30-25)/(18)*t;
    num=num+1;
end
for t = 18.5:0.5:305.5
    
    for r =1:5
        if t>t_region_vm(r,1)&t<t_region_vm(r,2)
            if r==1
                
                P1_x(num) = t;
                P1_y(num) = T_2(r)-b(r).*exp(-betaall(r,2)*(t-t_region_vm(r,1)));
                num=num+1;
            elseif r>1&&r<=5
                ychu(r) = y0(r-1)+k_interval_vm(r-1)*(t_region_vm(r,1)-t_region_vm(r-1,2));
                b(r)=T_2(r)-ychu(r);
                P1_x(num) = t;
                P1_y(num) = T_2(r)-b(r).*exp(-betaall(r,2)*(t-t_region_vm(r,1)));
                num=num+1;
           
            end
        elseif t>t_region_vm(r,2)&t<t_region_vm(r+1,1)
            k_interval_vm(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region_vm(r,2)-t_region_vm(r,1)));
            P1_x(num) = t;
            y0(r) = T_2(r)-b(r).*exp(-betaall(r,2)*(t_region_vm(r,2)-t_region_vm(r,1)));
            P1_y(num) = y0(r)+k_interval_vm(r)*(t-t_region_vm(r,2));
            num=num+1;
        end
    end
    
end
P1_x, P1_y

plot(P1_x,P1_y);

