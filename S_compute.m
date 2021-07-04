function [S,T_3,v,b,t_region,T_region,k_interval]=S_compute(f,betaall)
%%%%%%%%%%%%%%%%%%%%%给定一条染色体，生成五元组，检测并计算面积%%%%%%%%%%%%%%%%%%%%%%%
    s_region=[25,197.5;202.5,233;238,268.5;273.5,339.5;344.5,410.5];
	mid = [175;195;235;255];
	T_3 = mid+f(1:4,1)*10;
    T_3 = [T_3;25];
	v = (f(5,1)*17.5+82.5)/60;
	flag = 1;
	t_region=s_region/v;
    b(1)=-30+T_3(1);
	for r = 1:5
	    if r==1            
	        T_region(r,1)=30;
	        T_region(r,2)=b(1)*(1-exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1))))+T_region(r,1);
	        k_interval(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)));
	        T_region(r+1,1) = T_region(r,2) + k_interval(r)*(t_region(r+1,1)-t_region(r,2));
	        b(r+1) = T_3(r+1)-T_region(r+1,1);
	    elseif r<5
	        T_region(r,2)=b(r)*(1-exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1))))+T_region(r,1);
	        k_interval(r) = b(r)*betaall(r,2).*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)));
	        T_region(r+1,1) = T_region(r,2) + k_interval(r)*(t_region(r+1,1)-t_region(r,2));
	        b(r+1) = T_3(r+1)-T_region(r+1,1);
	    elseif r==5
	        T_region(r,2)=b(r)*(1-exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1))))+T_region(r,1);
	    end
	end

	%%%%%%%%%%%%%%判断斜率%%%%%%%%%%%%%%%%%%
	for r=1:5
		k_max(r)=abs(b(r)*betaall(r,2));
		if k_max(r)>3
			flag = 0;
		end
	end

	%%%%%%%%%%%%%判断150-190时间%%%%%%%%%%%%
	for r=1:4
		if 150>T_region(r,1)&150<T_region(r,2)
			t_150_190(1)=-1/betaall(r,2)*log((T_3(r)-150)/b(r))+t_region(r,1);
		
		elseif 150>T_region(r,2)&150<T_region(r+1,1)			
			t_150_190(1)=(150-T_region(r,2))/k_interval(r)+t_region(r,2);
		end	
		
		if 190>T_region(r,1)&190<T_region(r,2)
			t_150_190(2)=-1/betaall(r,2)*log((T_3(r)-190)/b(r))+t_region(r,1);
		
		elseif 190>T_region(r,2)&190<T_region(r+1,1)			
			t_150_190(2)=(190-T_region(r,2))/k_interval(r)+t_region(r,2);
		end	
	end
	if ~(60<t_150_190(2)-t_150_190(1)&120>t_150_190(2)-t_150_190(1))
		flag = 0;
	end
	%%%%%%%%%%%%%判断217时间%%%%%%%%%%%%%
	for r=1:4
		if 217>T_region(r,1)&217<T_region(r,2)
			t_217(1)=-1/betaall(r,2)*log((T_3(r)-217)/b(r))+t_region(r,1);							
	    %%%%%line%%%%%%%
		elseif 217>T_region(r,2)&217<T_region(r+1,1)			
			t_217(1)=(217-T_region(r,2))/k_interval(r)+t_region(r,2);
		end	
	end
	if 217>T_region(5,2)&217<T_region(5,1)
		t_217(2)=-1/betaall(5,2)*log((T_3(5)-217)/b(5))+t_region(5,1);
	end
	if ~(t_217(2)-t_217(1)>40&t_217(2)-t_217(1)<90)
		flag = 0;
	end
	%%%%%%%%%%%%判断峰值%%%%%%%%%%%%%%%%%%
	T_max=T_region(5,1);
	if ~(T_max>240&T_max<250)
		flag=0;
	end

	if flag==1
	%%%%%%%%%%%检测通过,计算面积%%%%%%%%%%%%%%
		for r=3:4
		    Si(r)=(T_region(r+1,1)+T_region(r,2))*(t_region(r+1,1)-t_region(r,2))*0.5;
		end
		for r=4:4
		    Sd(r) = T_3(r)*(t_region(r,2)-t_region(r,1))+(b(r)/betaall(r,2))*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)))-(b(r)/betaall(r,2));
		end
		for r=1:4
		    if 217>T_region(r,1)&217<T_region(r,2)

		        if r==3
		            S=T_3(r)*(t_region(r,2)-t_region(r,1))+(b(r)/betaall(r,2))*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)))-T_3(r)*(t_217(1)-t_region(r,1))-(b(r)/betaall(r,2))*exp(-betaall(r,2)*(t_217(1)-t_region(r,1)));
                    if S<0
                        S;
                    end
		            S=S+Si(3)+Sd(4)+Si(4);
		        elseif r==4
		            S=T_3(r)*(t_region(r,2)-t_region(r,1))+(b(r)/betaall(r,2))*exp(-betaall(r,2)*(t_region(r,2)-t_region(r,1)))-T_3(r)*(t_217(1)-t_region(r,1))-(b(r)/betaall(r,2))*exp(-betaall(r,2)*(t_217(1)-t_region(r,1)));
		            S=S+Si(4);
                    if S<0
                        S;
                    end
		        end
		    %%%%%%%%%%%%%line%%%%%%%%%%%%%%%%
		    elseif 217>T_region(r,2)&217<T_region(r+1,1)			
		        
		        if r==3
		            S=(T_region(r+1,1)+217)*(t_region(r+1,1)-t_217(1))*0.5;
		            S=S+Sd(4)+Si(4);
		        elseif r==4
		            S=(T_region(r+1,1)+217)*(t_region(r+1,1)-t_217(1))*0.5;

		        end
		    end	
        end
        if S<0
            S;
        end
		S=S-217*(t_region(5,1)-t_217(1)); 
        if S<0
            S;
        end
		for r = 5:5
            S_right=T_3(r)*(t_217(2)-t_region(r,1))+(b(r)/betaall(r,2))*exp(-betaall(r,2)*(t_217(2)-t_region(r,1)))-(b(r)/betaall(r,2))-217*(t_217(2)-t_region(r,1));            
        end

 
        
    else
        S=inf;
        S_right=inf;
    end