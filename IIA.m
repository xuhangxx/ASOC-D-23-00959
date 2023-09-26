function []=IIA()


for iter=1:11
    
    clearvars -except iter;
load 2000input.mat;
global history_list_1;
global num_history_list_1;
load one2one.mat;
load history.mat;
for i=1:n
    COEC(3,i)=deg2rad(COEC(3,i));
    COEC(4,i)=deg2rad(COEC(4,i));
    COEC(6,i)=deg2rad(COEC(6,i));
end
for i=1:m
    COET(3,i)=deg2rad(COET(3,i));
    COET(4,i)=deg2rad(COET(4,i));
    COET(6,i)=deg2rad(COET(6,i));
end
tic;
final=cell(n,1);                                      
N=1;
sum_fuel=zeros(1,N);
min_fuel=1e5*ones(1,N);
VOC0=VOC;
COEC0=COEC;
RTOM0=RTOM;
COET0=COET;
global num_cost;
num_cost=0;
num_iter=1;

for num_iter=1:N
    if(num_iter<=(N+1))
        AA_List=[];
        for j=1:n
            AA_List=[AA_List,j];
        end
    else
        AA_List=randperm(n);
    end
    for i=1:n
        COEC(:,i)=COEC0(:,AA_List(i));
        RTOM(:,i)=RTOM0(:,AA_List(i));
        VOC(:,i)=VOC0(:,AA_List(i));
        COET=COET0;
    end
    kdv=cal_k(n,m,TOT);
    Mu=3.98600441500e14;    
    Final=cell(n,1);
    outlist=zeros(m,1);
    %one2one=cell(m,n);
    %[list_serving_codeold,list_serving_code,one2one,num_serving_T,num_fenpei]=onetoone(n,m,RTOM,VOC,TOT,COEC,COET);  
    num_cost=num_cost+48;
    toc;
    ttime=toc-tic;
    %load one2one09112.mat;
    list_dis=1e10*ones(1,m);                      
    ifcho_m=ones(1,m);                            
    numofc=zeros(1,n);                             
    if(num_iter==1)
        [list_dis]=cal_meandis(COEC,one2one,COET,ifcho_m,n,m,kdv,TOT,VOC);  
    else
        [list_dis]=cal_meandis2(COEC,COET,ifcho_m,n,m,kdv,TOT,VOC);  
    end
    [list_cho]=cal_cho(list_dis,n,m,ifcho_m);       
    VOC0=VOC;
    delta=mean(VOC)/m*n/10;
    list_his=-1*ones(n,m);         
    num_m_ben=size(list_cho,2);       
    list_fuel=zeros(num_m_ben,n);
    for i=1:num_m_ben         
        for j=1:n
            if(one2one{list_cho(i),j}==-1)
                list_fuel(i,j)=-1;
            else
                list_fuel(i,j)=one2one{list_cho(i),j}(3,1);
            end
        end
    end
    %[list_fuel,finalx]=cal_cost(list_cho,COEC,COET,n,num_m_ben,TOT,list_his,VOC,VOC0,numofc,RTOM);     
    list_cost=list_fuel;                               
    num_pai=0;                                         
    num_m_got=0;                                      
    num_m_cal = m;                                    
    flag_liupei=1;
    COEC_2=COEC;
    %list_jinji=ones(n,m);
    while(num_m_cal>0 && flag_liupei)
        tarofc=-1*ones(1,n);                             
        num_pai=num_pai+1;
        flag_liupei=0;
        list_m_cho=zeros(num_m_ben,1);
        for c=1:n              
            flagc=1;
            min_cost=1e10;
            min_number=0;
            m_youxiao=0;
            for i=1:num_m_ben
                if(list_cost(i,c)>-0.1 && list_serving_code(list_cho(i),c)==1)
                    m_youxiao= m_youxiao+1;
                    if(list_cost(i,c)<min_cost)
                        min_cost=list_cost(i,c);
                        %m_youxiao= m_youxiao+1;
                        min_number=i;   
                    end
                end
            end
            if(min_cost>1e9)
                flagc=0;
                continue;
            end
            if(list_fuel(min_number,c)<VOC(c)&&list_fuel(min_number,c)>0)
                if(outlist(list_cho(min_number))~=0)
                    bef_number=outlist(list_cho(min_number));             
                    outlist(list_cho(min_number))=c;
                    tarofc(c)=list_cho(min_number);
                    numofc(c)=numofc(c)+1;
                    tarofc(bef_number)=-1;
                    numofc(bef_number)=numofc(bef_number)-1;
                    if(num_m_ben>1&& m_youxiao>1)%
                        min_cost2=1e10;
                        for k=1:num_m_ben
                            if(list_serving_code(list_cho(k),c)==1 && list_cost(k,c)<min_cost2 && list_cost(k,c)>0 && k~=min_number && list_cost(k,c)< VOC(c)+list_fuel(k,c))
                                min_cost2=list_cost(k,c);
                            end
                        end
                        for j=1:n
                            if(list_cost(min_number,j)>0)
                                list_cost(min_number,j)=list_cost(min_number,j)+min_cost2-min_cost+delta;      
                            end
                        end
                    elseif(m_youxiao>0)
                        for j=1:n
                            if(list_cost(min_number,j)>0)
                                list_cost(min_number,j)=list_cost(min_number,j)+delta;     
                            end
                        end
                    end
                    VOC(c)=VOC(c)-list_fuel(min_number,c);
                    flag_recho=0;
                    cost_2=1e10;                                            
                    number_2=-1;
                    m_youxiao=0;
                    for u=1:num_m_ben
                        if(list_serving_code(list_cho(u),bef_number)==1 && list_cost(u,bef_number)>-0.1 && list_cost(u,bef_number)< VOC(bef_number)+list_fuel(u,bef_number) && list_cost(u,bef_number)< VOC(bef_number)+list_fuel(u,bef_number))
                            m_youxiao= m_youxiao+1;
                            if(list_cost(u,bef_number)<cost_2)
                                flag_recho=1;
                                cost_2=list_cost(u,bef_number);
                                number_2=u;
                            end
                        end
                    end
                    if(flag_recho==1 &&list_fuel(number_2,bef_number)<VOC(bef_number)+list_fuel(min_number,bef_number))                      
                        if(outlist(list_cho(number_2))~=0 )                                        
                            while (outlist(list_cho(number_2))~=0 )
                                if(~flag_recho)
                                    break;
                                end
                                bef_number_2=outlist(list_cho(number_2));
                                
                                outlist(list_cho(number_2))=bef_number;       
                                tarofc(bef_number)=list_cho(number_2);
                                numofc(bef_number)=numofc(bef_number)+1;
                                numofc(bef_number_2)=numofc(bef_number_2)-1;
                                tarofc(bef_number_2)=-1;
                                if(num_m_ben>1&& m_youxiao>1)
                                    min_cost2=1e10;
                                    for k=1:num_m_ben
                                        if(list_serving_code(list_cho(k),bef_number)==1 && list_cost(k,bef_number)<min_cost2 && list_cost(k,bef_number)>0 && k~=number_2 && list_cost(k,bef_number)< VOC(bef_number)+list_fuel(k,bef_number))
                                            min_cost2=list_cost(k,bef_number);
                                        end
                                    end
                                    for j=1:n
                                        if(list_cost(number_2,j)>0)
                                            list_cost(number_2,j)=list_cost(number_2,j)+min_cost2-cost_2+delta;      
                                        end
                                    end
                                elseif(m_youxiao>0)
                                    for j=1:n
                                        if(list_cost(number_2,j)>0)
                                            list_cost(number_2,j)=list_cost(number_2,j)+delta;     
                                        end
                                    end
                                end
                                VOC(bef_number)=VOC(bef_number)+list_fuel(min_number,bef_number)-list_fuel(number_2,bef_number);
                                
                                bef_number=bef_number_2;
                                min_number=number_2;
                               
                                flag_recho=0;
                                cost_2=1e10;                                  
                                m_youxiao=0;
                                for u=1:num_m_ben
                                    if( list_serving_code(list_cho(u),bef_number_2)==1&& list_cost(u,bef_number_2)>-0.1 && list_cost(u,bef_number_2)< VOC(bef_number_2)+list_fuel(u,bef_number_2) && list_cost(u,bef_number_2)< VOC(bef_number_2)+list_fuel(u,bef_number_2))
                                        m_youxiao= m_youxiao+1;
                                        if(list_cost(u,bef_number_2)<cost_2)
                                            flag_recho=1;
                                            if(879<cost_2&&cost_2<880)
                                                a=1;
                                            end
                                            cost_2=list_cost(u,bef_number_2);
                                            number_2=u;
                                        end
                                    end
                                end
                                if(flag_recho==1 && list_fuel(number_2,bef_number)<VOC(bef_number)+list_fuel(min_number,bef_number))                        
                                    if(outlist(list_cho(number_2))==0 )                                         
                                      
                                        outlist(list_cho(number_2))=bef_number;
                                        list_m_cho(number_2)=1;
                                        tarofc(bef_number)=list_cho(number_2);
                                        numofc(bef_number)=numofc(bef_number)+1;
                                        if(num_m_ben>1&& m_youxiao>1)
                                            min_cost2=1e10;
                                            for k=1:num_m_ben
                                                if(list_serving_code(list_cho(k),bef_number)==1 && list_cost(k,bef_number)<min_cost2 && list_cost(k,bef_number)>0 && k~=number_2 && list_cost(k,bef_number)< VOC(bef_number)+list_fuel(k,bef_number))
                                                    min_cost2=list_cost(k,bef_number);
                                                end
                                            end
                                            for j=1:n
                                                if(list_cost(number_2,j)>0)
                                                    list_cost(number_2,j)=list_cost(number_2,j)+min_cost2-cost_2+delta;      
                                                end
                                            end
                                        elseif(m_youxiao>0)
                                            for j=1:n
                                                if(list_cost(number_2,j)>0)
                                                    list_cost(number_2,j)=list_cost(number_2,j)+delta;      %
                                                end
                                            end
                                        end
                                        VOC(bef_number)=VOC(bef_number)+list_fuel(min_number,bef_number)-list_fuel(number_2,bef_number);
                                        num_m_got=num_m_got+1;
                                        break;
                                    end
                                else
                                    VOC(bef_number)=VOC(bef_number)+list_fuel(min_number,bef_number);
                                    break;
                                end
                            end
                        else
                            outlist(list_cho(number_2))=bef_number;
                            list_m_cho(number_2)=1;
                            tarofc(bef_number)=list_cho(number_2);
                            numofc(bef_number)=numofc(bef_number)+1;
                            if(num_m_ben>1&& m_youxiao>1)
                                min_cost2=1e10;
                                for k=1:num_m_ben
                                    if(list_serving_code(list_cho(k),bef_number)==1 && list_cost(k,bef_number)<min_cost2 && list_cost(k,bef_number)>0 && k~=number_2 && list_cost(k,bef_number)< VOC(bef_number)+list_fuel(k,bef_number))
                                        min_cost2=list_cost(k,bef_number);
                                    end
                                end
                                for j=1:n
                                    if(list_cost(number_2,j)>0)
                                        list_cost(number_2,j)=list_cost(number_2,j)+min_cost2-cost_2+delta;      
                                    end
                                end
                            elseif(m_youxiao>0)
                                for j=1:n
                                    if(list_cost(number_2,j)>0)
                                        list_cost(number_2,j)=list_cost(number_2,j)+delta;      
                                    end
                                end
                            end
                            VOC(bef_number)=VOC(bef_number)+list_fuel(min_number,bef_number)-list_fuel(number_2,bef_number);
                            num_m_got=num_m_got+1;
                        end
                    else
                        VOC(bef_number)=VOC(bef_number)+list_fuel(min_number,bef_number);
                    end
                else
                    outlist(list_cho(min_number))=c;
                    list_m_cho(min_number)=1;
                    tarofc(c)=list_cho(min_number);
                    if(num_m_ben>1&& m_youxiao>1)
                        min_cost2=1e10;
                        for k=1:num_m_ben
                            if(list_serving_code(list_cho(k),c)==1 && list_cost(k,c)<min_cost2 && list_cost(k,c)>0 && k~=min_number && list_cost(k,c)< VOC(c)+list_fuel(k,c))
                                min_cost2=list_cost(k,c);
                            end
                        end
                        for j=1:n
                            if(list_cost(min_number,j)>0)
                                list_cost(min_number,j)=list_cost(min_number,j)+min_cost2-min_cost+delta;     
                            end
                        end
                    elseif(m_youxiao>0)
                        for j=1:n
                            if(list_cost(min_number,j)>0)
                                list_cost(min_number,j)=list_cost(min_number,j)+delta;     
                            end
                        end
                    end
                    numofc(c)=numofc(c)+1;
                    VOC(c)=VOC(c)-list_fuel(min_number,c);
                    num_m_got=num_m_got+1;
                end
            else
                flagc=0;
            end
        end
        for i=1:n
            if(tarofc(i)~=-1)
                COEC_2(:,i)=COET(:,tarofc(i));
                list_his(i,num_pai)=tarofc(i);
            end
        end
        for i=1:num_m_ben
            if(list_m_cho(i)==1)
                ifcho_m(list_cho(i))=0;
                flag_liupei=1;
            end
        end
        num_m_cal=num_m_cal-sum(list_m_cho);
        if( num_m_cal>0 && flag_liupei)
            if(num_iter==1)
                [list_dis]=cal_meandis(COEC_2,one2one,COET,ifcho_m,n,m,kdv,TOT,VOC); 
            else
                [list_dis]=cal_meandis2(COEC_2,COET,ifcho_m,n,m,kdv,TOT,VOC);  
            end
            [list_cho]=cal_cho(list_dis,n,num_m_cal,ifcho_m);              
            delta=mean(VOC)/m*n/10;
            num_m_ben=size(list_cho,2);
            [list_fuel,finalx]=cal_cost(list_cho,COEC,COET,n,num_m_ben,TOT,list_his,VOC,VOC0,numofc,RTOM);     
            list_cost=list_fuel;
        end
    end
    [final]=cal_cost2(COEC,COET,n,m,TOT,list_his,VOC0,numofc,RTOM,AA_List);     
    sum_fuel(num_iter)=0;
    for i=1:n
        for j=1:size(final{i,1},2)
            sum_fuel(num_iter)=sum_fuel(num_iter)+final{i,1}(3,j);%Final{i,1}(10,j);
        end
    end
    if(num_iter==1)
        min_fuel(num_iter)=sum_fuel(num_iter);
    else
        if(sum_fuel(num_iter)<min_fuel(num_iter-1))
            min_fuel(num_iter)=sum_fuel(num_iter);
        else
            min_fuel(num_iter)=min_fuel(num_iter-1);
        end
    end
    toc;
    ttime=toc-tic;
end
toc;
ttime=toc-tic;
tic;
VOC=VOC0;
COEC=COEC0;
RTOM=RTOM0;
COET=COET0;
%% parameter
com_parameter.threshold=0.0052;                             
com_parameter.KT=100;
com_parameter.L=100;
com_parameter.K=0.9;
com_parameter.KT2=10;
com_parameter.L2=10;
com_parameter.K2=0.5;
Mu=3.98600441500e14;    %miu,m^3/s^2
final_0=final;

num_m2=0;
flag_fenpei=1;         
TT=100;
T=100;
number_iter=1;
FLAG1=1;
num_iter_1=log(0.001/T)/log(com_parameter.K);        
%5s
len(1)=0;
one_hist=cell(1,1);
two_hist=zeros(1,1);
three_hist=zeros(1,1);
std_hist=zeros(1,m);
final_pre=final_0;
fzq=zeros(com_parameter.L,m);
city=zeros(1,m);
for i=1:n
    for j=1:size(final{i,1},2)
        city(final{i,1}(1,j))=i;
    end
end
city0_1=[];
num_0city=0;                 %零变量个数
for j=1:m
    if(city(j)==0)
        num_0city=num_0city+1;
        city0_1(num_0city)=j;  %零变量目标编号
    end
end
m_ran=0;                     %零变量个数
if(num_0city<(m-1))
    m_ran=m-num_0city-1;
else
    m_ran=m-num_0city;
end
flag_new_1=0;
num_pro=2;
num_begin_now=1;

for i=num_begin_now:(com_parameter.L-1)
    city0=city0_1;
    city_copy=city;
    numchange=m-num_0city;
    change_sum=floor(numchange*rand()+1);
    for j=1:change_sum
        m_cho=floor(m*rand()+1);
        while(city_copy(m_cho)==0)
            m_cho=floor(m*rand()+1);
        end
        city_copy(m_cho)=0;
        city0(num_0city+j)=m_cho;  %零变量目标编号
    end
    change_sum=floor((change_sum+num_0city)*rand()+1);
    j_list=city0(randperm(numel(city0),change_sum));
    flag_no=0;
    for k=1:size(j_list,2)
        if(sum(list_serving_code(j_list(k),:))==0)
            flag_no=1;
        end
    end
    while (flag_no)
        j_list=city0(randperm(numel(city0),change_sum));
        for k=1:size(j_list,2)
            if(sum(list_serving_code(j_list(k),:))==0)
                flag_no=1;
            end
        end
    end
    for k=1:size(j_list,2)
        if(city_copy(j_list(k))==0)
            list_serving_code_t=list_serving_code(j_list(k),:);
            flag_gotn=1;
            while(flag_gotn)
                n_cho=floor(n*rand()+1);
                if(list_serving_code_t(1,n_cho)==1)
                    city_copy(j_list(k))=n_cho;
                    flag_gotn=0;
                end
            end
        end
    end
    fzq(i,:)=city_copy;
end
num_begin_now=(com_parameter.L-1);
global max_mgot;
global number_ini;
global number_max;
number_ini=0;
number_max=0;
for i=num_begin_now:com_parameter.L
    fzq(i,:)=city;
end
max_mgot=m-num_0city+num_pro-2;
number_ini=num_cost;
number_m_inigot=max_mgot;
num1=zeros(1,1);
rate_deltae=100*T*com_parameter.K.^(round(num_iter_1/2));
d_len1=1;
fzq_best=city;
cout_dead=0;
toc;
ttime=toc-tic;
four_hist=zeros(1,1);
five_hist=zeros(1,1);
load one2one0.mat;
while (T > 0.001|| cout_dead<50)
    one_hist{number_iter,1}=fzq;
    for i=1:m
        std_hist(number_iter,i)=std(one_hist{number_iter,1}(:,i));
    end
    std_hist(number_iter,m+1)=mean(std_hist(number_iter,1:m));
    if(d_len1>=0)
        [nzq,len(number_iter),fzq_best,num_m2(number_iter),two_hist(number_iter),three_hist(number_iter)]=got_pinggu(fzq,m,Mu,n,TOT,RTOM,VOC,COEC,COET,final_pre,com_parameter);
        if number_iter>1
            d_len1=(len(number_iter)-len(number_iter-1))/len(number_iter-1);
            if d_len1<0
                d_len1=1e4;
            end
        else
            d_len1=1e4;
        end
        [fzq,four_hist(number_iter),five_hist(number_iter)]=generate_newzq(nzq,fzq_best,m,list_serving_codeold,list_serving_code,num_serving_T,rate_deltae,T,Mu,n,TOT,RTOM,VOC,COEC,COET,final_pre,com_parameter,num_m2(number_iter));
    end
    TT(number_iter)=T;
    T=T*com_parameter.K;
    if(number_iter>1)
        if(len(number_iter-1)==len(number_iter))
            cout_dead=cout_dead+1;
        else
            cout_dead=0;
        end
    end
    number_iter=number_iter+1;
    jindu2=floor((number_iter-1)/num_iter_1*100);
    %disp(jindu2);
    if(jindu2>98)
        a=1;
    end
    final_pre=final_0;
    disp(T);
    toc;
    ttime=toc-tic;
end
city=fzq_best;
[~,~,~,~,final,~]=func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter);

toc;
ttime=toc-tic;
eval(['save(''200iia_result',num2str(iter),'.mat'');'])
end


flag=1;

function [nzq,maxlen,fzq_best,num_mlast,two_hist,three_hist]=got_pinggu(fzq,m,Mu,n,TOT,RTOM,VOC,COEC,COET,final_pre,com_parameter)

global max_mgot;
global num_cost;
global number_max;
nzq(1,:)=fzq(1,:);
for i=1:size(fzq,1)
    city=fzq(i,:);
    [lenzq(i),~,~,flagzq(i),final,num_m(i)]=func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter);
    if(flagzq(i)==1 && num_m(i)>max_mgot)
        max_mgot=num_m(i);
        number_max=num_cost;
        disp('max:');
        disp(max_mgot);
        disp('number_max:');
        disp(num_cost);
    end
end
two_hist=0;
maxlen=0;
minlen=1e10;
num_youxiao=0;
num_max=0;
for i=1:size(fzq,1)
    if(flagzq(i)~=0)
        num_youxiao=num_youxiao+1;
        two_hist=two_hist+lenzq(i);
        if(maxlen<lenzq(i))
            maxlen=lenzq(i);
        end
        if(minlen>lenzq(i))
            minlen=lenzq(i);
        end
        if(num_m(i)>num_max)
             num_max=num_m(i);
        end
    end
end
number_minlen=find(lenzq==maxlen);
two_hist=two_hist/num_youxiao;
three_hist=num_youxiao;
if(size(number_minlen,2)>1)
    fzq_best=fzq(number_minlen(1),:);   
    num_mlast=num_m(number_minlen(1));
else
    fzq_best=fzq(number_minlen,:);  
    num_mlast=num_m(number_minlen);
end
if(size(number_minlen,2)==0)
    num_mlast=0;
    fzq_best=zeros(1,m);
end
if(num_max~=num_mlast)
    disp('eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee');
end
%%  update  one
for i=1:length(lenzq)
    if(flagzq(i)~=0)
        if(num_youxiao<com_parameter.L)
            dddd=1;
        end
        if(num_youxiao<com_parameter.L/2)
            fitness(i,1)=1;
        else
            rate=(num_youxiao-com_parameter.L/2)*2/com_parameter.L;
            fitness(i,1)=(1-rate*(maxlen-lenzq(i))/(maxlen-minlen+0.001));
        end
    end
end


nn=0;
if(maxlen==minlen)
    for i=1:size(fzq,1)
        if(flagzq(i)~=0)
            if(rand()>0.5)
                nn=nn+1;
                nzq(nn,:)=fzq(i,:);
            end
        end
    end
else
    for i=1:size(fzq,1)
        if(flagzq(i)~=0)
            if(fitness(i,1)>=rand())
                nn=nn+1;
                nzq(nn,:)=fzq(i,:);
            end
        end
    end
end
if(flagzq(1)==0&& size(nzq,1)==1)
    nzq(1,:)=zeros(1,m);
end
if(minlen>1e9)
    minlen=-1;
end
return

%generate new population (GA)
function [fzq,four_hist,five_hist]=generate_newzq(nzq,fzq_best,m,list_serving_codeold,list_serving_code,num_serving_T,rate_deltae,T,Mu,n,TOT,RTOM,VOC,COEC,COET,final_pre,com_parameter,num_m)
num_zq=size(nzq,1);
five_hist=num_zq;
num_gen=0;
if(num_zq>1)
    num_add=0;
    for i=1:com_parameter.L
        [baby1(i,:),baby2(i,:),father1(i,:),father2(i,:)]=generate_newbaby(num_zq,nzq,m,list_serving_codeold,num_serving_T,num_m);
        city=father1(i,:);
        [len_f1(i,1),sumt,sumfuel,flag_f1(i,1),final,num_mi]=func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter);
        city=baby1(i,:);
        [len_b1(i,1),sumt,sumfuel,flag_b1(i,1),final,num_mi]=func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter);
        city=father2(i,:);
        [len_f2(i,1),sumt,sumfuel,flag_f2(i,1),final,num_mi]=func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter);
        city=baby2(i,:);
        [len_b2(i,1),sumt,sumfuel,flag_b2(i,1),final,num_mi]=func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter);
        delta_e=rate_deltae*(len_f1(i,1)-len_b1(i,1))/len_f1(i,1);
        if (flag_f1(i,1)==1 && flag_b1(i,1)==1)
            if delta_e<0
                nzq=[nzq;baby1(i,:)];
                len_b1(i,1)=-2;
                num_add=num_add+1;
                if(num_add==com_parameter.L-num_zq)
                    break;
                end
            else
                %% update two
                if (exp(-delta_e/T)>rand())
                    nzq=[nzq;baby1(i,:)];
                    len_b1(i,1)=-2;
                    num_add=num_add+1;
                    if(num_add==com_parameter.L-num_zq)
                        break;
                    end
                end
            end
        elseif(flag_f1(i,1)==0)
            nzq=[nzq;baby1(i,:)];
            len_b1(i,1)=-2;
            num_add=num_add+1;
            if(num_add==com_parameter.L-num_zq)
                break;
            end
        end
        delta_e=rate_deltae*(len_f2(i,1)-len_b2(i,1))/len_f2(i,1);
        if (flag_f2(i,1)==1 && flag_b2(i,1)==1)
            if delta_e<0
                nzq=[nzq;baby2(i,:)];
                len_b2(i,1)=-2;
                num_add=num_add+1;
                if(num_add==com_parameter.L-num_zq)
                    break;
                end
            else
                if (exp(-delta_e/T)>rand())
                    nzq=[nzq;baby2(i,:)];
                    len_b2(i,1)=-2;
                    num_add=num_add+1;
                    if(num_add==com_parameter.L-num_zq)
                        break;
                    end
                end
            end
        elseif(flag_f2(i,1)==0)
            nzq=[nzq;baby2(i,:)];
            len_b2(i,1)=-2;
            num_add=num_add+1;
            if(num_add==com_parameter.L-num_zq)
                break;
            end
        end
    end
    num_zq=size(nzq,1);
    while(num_zq<com_parameter.L)
        xuhao_coft=floor(rand*num_zq+1);
        nzq=[nzq;nzq(xuhao_coft,:)];
        num_zq=size(nzq,1);
    end
elseif(num_zq>0)
    a=1;
    for i=2:com_parameter.L
        nzq(i,:)=nzq(1,:);
    end
else
    AA=rand(com_parameter.L,m);
    for i=1:com_parameter.L
        for j=1:m
            num_oft=num_serving_T(j);
            xuhao_coft=floor(AA(i,j)*num_oft+1);
            list_serving_code_t=list_serving_codeold{j,1};
            city(j)=list_serving_code_t(1,xuhao_coft);
        end
        nzq(i,:)=city;
    end
end
four_hist=num_gen;
%% update three
if(num_zq>1)
    for i=round(num_zq/2):num_zq
        num_0city=0;
        for j=1:m
            if(nzq(i,j)==0)
                num_0city=num_0city+1;
                city0(num_0city)=j;
            end
        end
        aa=0;
        for j=1:num_0city
            aa=aa+sum(list_serving_code(city0(j),:));
        end
        if(num_0city>0&& aa>0)
            j=city0(floor(num_0city*rand()+1));
            while (sum(list_serving_code(j,:))==0)
                j=city0(floor(num_0city*rand()+1));
            end
            if(nzq(i,j)==0)
                list_serving_code_t=list_serving_code(j,:);
                flag_gotn=1;
                while(flag_gotn)
                    n_cho=floor(n*rand()+1);
                    if(list_serving_code_t(1,n_cho)==1)
                        nzq(i,j)=n_cho;
                        flag_gotn=0;
                    end
                end
            end
        end
    end
end
if(num_zq>com_parameter.L)
    nzq=nzq(num_zq-com_parameter.L+1:num_zq,:);  
end
fzq=nzq;
if(~isempty(fzq_best))
    fzq(1,:)=fzq_best; 
end
return

function [baby1,baby2,father1,father2]=generate_newbaby(num_fzq,nzq,m,list_serving_codeold,num_serving_T,num_m)
num_fzqcho=randperm(num_fzq);
baby1=nzq(num_fzqcho(1),:);
baby2=nzq(num_fzqcho(2),:);
father1=baby1;
father2=baby2;
%%%cross
W=ceil(m/10);  
p=floor(1+rand()*(m-W+1));%% [0-1）
for i=1:W
    x=find(baby1==baby2(1,p+i-1));   
    y=find(baby2==baby1(1,p+i-1));
    max_iter=100;
    while(isempty(x)||isempty(y))
        max_iter=max_iter-1;
        if(max_iter<0)
            break;
        end
        p=ceil(rand()*(m-W+1));
        x=find(baby1==baby2(1,p+i-1));
        y=find(baby2==baby1(1,p+i-1));
    end
    if(max_iter<0)
        break;
    end
    temp=baby1(1,p+i-1);  
    baby1(1,p+i-1)=baby2(1,p+i-1);
    baby2(1,p+i-1)=temp;
    x2=x(floor(1+size(x,2)*rand()));
    temp=baby1(1,x2);
    y2=y(floor(1+size(y,2)*rand()));
    baby1(1,x2)=baby2(1,y2);
    baby2(1,y2)=temp;
end
%%mutation
p1=floor(1+m*rand());
p2=floor(1+m*rand());
while(p1==p2 || num_serving_T(p1)<2 || num_serving_T(p2)<2)
    p1=floor(1+m*rand());
    p2=floor(1+m*rand());
end
tmp=baby1(p1);
tmp2=list_serving_codeold{p1,1}(1,floor(rand(1)*num_serving_T(p1)+1));
while tmp==tmp2
    tmp2=list_serving_codeold{p1,1}(1,floor(rand(1)*num_serving_T(p1)+1));
end
baby1(p1)=tmp2;
tmp=baby1(p2);
tmp2=list_serving_codeold{p2,1}(1,floor(rand(1)*num_serving_T(p2)+1));
while tmp==tmp2
    tmp2=list_serving_codeold{p2,1}(1,floor(rand(1)*num_serving_T(p2)+1));
end
baby1(p2)=tmp2;
tmp=baby2(p1);
tmp2=list_serving_codeold{p1,1}(1,floor(rand(1)*num_serving_T(p1)+1));
while tmp==tmp2
    tmp2=list_serving_codeold{p1,1}(1,floor(rand(1)*num_serving_T(p1)+1));
end
baby2(p1)=tmp2;
tmp=baby2(p2);
tmp2=list_serving_codeold{p2,1}(1,floor(rand(1)*num_serving_T(p2)+1));
while tmp==tmp2
    tmp2=list_serving_codeold{p2,1}(1,floor(rand(1)*num_serving_T(p2)+1));
end
baby2(p2)=tmp2;
return

%without perms，
function list_bh=bh2list(tarlist_ofc,num_t_ofc,bh_oflist)
if(num_t_ofc == 1)
    list_bh = tarlist_ofc(1);
end
if(num_t_ofc == 0)
    list_bh = [];
end
if(num_t_ofc > 1)
    list_bh = zeros(1,num_t_ofc);
    for i = 1:num_t_ofc
        umap(1,i) =i;
        umap(2,i) = tarlist_ofc(i);
        val(i) = i;
    end
    factor(1) = 1;
    for i = 2:num_t_ofc-1
        factor(i) = factor(i-1) *i;
    end
    bh_oflist = bh_oflist-1;
    for i = 1:num_t_ofc
        index = floor(bh_oflist / factor(num_t_ofc - i)) + i;
        temp = val(index);
        while(index>i)
            val(index) = val(index - 1);
            index = index-1;
        end
        val(i) = temp;
        bh_oflist = bh_oflist - floor(bh_oflist/factor(num_t_ofc - i)) * factor(num_t_ofc - i);
        if (bh_oflist == 0)
            break;
        end
    end
    for i = 1:num_t_ofc
        list_bh(i) = umap(2,val(i));
    end
end
return


function [list_dis]=cal_meandis(COEC_fk,one2one,COET_fk,ifcho_m,n,m,kdv,TOT,VOC)   
list_dis=1e10*ones(1,m);
list_cho=zeros(1,m);
ifcho_m2=ifcho_m;
for j=1:n
    tar_j=-1;
    cost_min=inf;
    rvc_c=COE2RV(COEC_fk(:,j));
    for i=1:m
        if(ifcho_m2(i))
            
            rvc_t=COE2RV(COET_fk(:,i));
            hc = cross(rvc_c(1:3), rvc_c(4:6));
            ht = cross(rvc_t(1:3), rvc_t(4:6));
            numda=COEC_fk(1)/COET_fk(1,i);
            thetaym=acos_2(dot(hc,ht)/(norm(hc)*norm(ht)));
            [theta1,theta2]=get2points(COEC_fk,COET_fk(:,i));
            thetaph=theta2-theta1;
            if(thetaph<-pi)
                thetaph=thetaph+2*pi;
            elseif(thetaph>pi)
                thetaph=thetaph-2*pi;
            end
            cost=2*abs(sqrt(2-(1+thetaph/2/pi).^(-2/3))-1)+2*sin(thetaym/2)+abs(sqrt(2/(numda*numda+numda))-sqrt(1/numda))+abs(1-sqrt(2*numda/(numda+1)));
            if(cost_min>cost)
                cost_min=cost;
                tar_j=i;
            end
        end
    end
    if(tar_j>0)
        list_dis(tar_j)=1;
        ifcho_m2(tar_j)=0;
    end
end
return

function [list_dis]=cal_meandis2(COEC_fk,COET_fk,ifcho_m,n,m,kdv,TOT,VOC)     
list_dis=1e30*ones(1,m);

list=randperm(sum(ifcho_m));
j=1;
for i=1:m
     if(ifcho_m(i))
         list_dis(i)=list(j);
         j=j+1;
     end
end
return

function [list_cho]=cal_cho(list_dis,n,m,ifcho_m)    
[A,ms]=sort(list_dis);
list_cho=[];
if (m<=n)
    num_xu=1;
    for i=1:m
        if(num_xu>m)
            break;
        end
        if(ifcho_m(ms(i))&& list_dis(ms(i))<1e9)
            num_xu=num_xu+1;
            list_cho=[list_cho,ms(i)];
        end
    end
else
    num_xu=1;
    for i=1:m
        if(num_xu>n)
            break;
        end
        if(ifcho_m(ms(i)) && list_dis(ms(i))<1e9)
            num_xu=num_xu+1;
            list_cho=[list_cho,ms(i)];
        end
    end
end

return

function [list_cost,finalx]=cal_cost(list_cho,COEC,COET,n,m,TOT,list_his,VOC,VOC0,numofc,RTOM)
global num_cost;
Mu=3.98600441500e14;    %miu,
list_cost=-1*ones(m,n);
finalx=cell(n,m);
for i=1:m
    for j=1:n
        if (numofc(j)==0)
            cTOT=TOT;
            cVOC=VOC0(j);
        else
            cTOT=TOT/(numofc(j)+1);
            cVOC=VOC0(j)/(numofc(j)+1);
        end
        coec=COEC(:,j);
        num_t_ofc=(numofc(j)+1);
        k_got=1;
        k_num=1;
        while(k_got<num_t_ofc)
            if(list_his(j,k_num)~=-1)
                tar=list_his(j,k_num); 
                city_2(k_got)=tar;
                k_got=k_got+1;
            end
            k_num=k_num+1;
        end
        city_2(k_got)=list_cho(i);
        [len,flagt,sumtc,finalx{i,j}]=func3(j,city_2,num_t_ofc,RTOM,cTOT,cVOC,coec,COET);
        num_cost=num_cost+1;
        if(flagt~=0)
            list_cost(i,j)=len+VOC(j)-VOC0(j);
        end
    end
end
return

function [final]=cal_cost2(COEC,COET,n,m,TOT,list_his,VOC0,numofc,RTOM,AA_List)
Mu=3.98600441500e14;    %miu,
final=cell(n,1);
for j=1:n
    if (numofc(j)==1)
        cTOT=TOT;
        cVOC=VOC0(j);
    else
        cTOT=TOT/(numofc(j));
        cVOC=VOC0(j)/(numofc(j));
    end
    coec=COEC(:,j);
    num_t_ofc=numofc(j);
    if(num_t_ofc==0)
        continue;
    end
    k_got=1;
    k_num=1;
    while(k_got<(num_t_ofc+1))
        if(list_his(j,k_num)~=-1)
            tar=list_his(j,k_num);   
            city_2(k_got)=tar;
            k_got=k_got+1;
        end
        k_num=k_num+1;
    end
    [len,flagt,sumtc,final{AA_List(j),1}]=func3(j,city_2,num_t_ofc,RTOM,cTOT,cVOC,coec,COET);
end
return


function a = Accel(x,xk)
GM_Earth = 398600.4415;
r_E = 6378.1363;                              
J2 = 1082.6355e-6;                           
r = sqrt((x(1) + xk(1))^2 + (x(2) + xk(2))^2 + (x(3) + xk(3))^2);

a(1) = -GM_Earth * (x(1) + xk(1))/ r^3 * (1 + 1.5 * J2 * (r_E / r)^2 * (1 - 5 * (x(3) + xk(3))^2 / r^2));
a(2) = a(1) * (x(2) + xk(2)) / (x(1) + xk(1));
a(3) = -GM_Earth * (x(3) + xk(3))/ r^3 * (1 + 1.5 * J2 * (r_E / r)^2 * (3 - 5 * (x(3) + xk(3))^2 / r^2));
return


function [len,sumt,sumfuel,flag,final,num_mi] = func1(m,city,RTOM,Mu,n,TOT,VOC,COEC,COET,final_pre,com_parameter)
global num_cost;
global history_list_1;
global num_history_list_1;
final_pre2=final_pre;
history_list=history_list_1;
flagi=1;
costfuel=0;                                             
costtime=0;                                             
VOC_after=VOC;
num_mi=0;
mode=2;
for c=1:n
    RVC=COE2RV(COEC(:,c));
    finalc=[];
    num=1;
    tarlist_ofc=find(city==c);
    if(isempty(tarlist_ofc))                   
        continue;
    end
    
    num_cost=num_cost+1;
    num_t_ofc=size(tarlist_ofc,2);
    num_shixu_c=1;
    for k=1:num_t_ofc
        num_shixu_c=num_shixu_c*k;
    end
    flag_history=0;
    if(~isempty(history_list{1,1}))                   
        for k=1:size(history_list,1)
            if(isequal(history_list{k,2},tarlist_ofc)&& history_list{k,1}==c)
                flag_history=1;
                costminc=history_list{k,3};
                if(costminc<VOC(c))
                    flagc=1;
                else
                    flagc=0;
                end
                costminc=history_list{k,3};
                sumtc=history_list{k,4};
                finalc=history_list{k,6};
                if(flagc==1)
                    costfuel=costfuel+costminc;
                    VOC_after(c)=VOC(c)-costminc;
                    if(sumtc>costtime)
                        costtime=sumtc;
                    end
                    final_pre2{c,1}=finalc;
                else
                    flagi=0;
                    break;
                end
            end
        end
    end
    if(flag_history==0)
        if(num_t_ofc>5)
            flagi=0;
            break;
        end
        num_history_list_1=num_history_list_1+1;
        history_list_1{num_history_list_1,1}=c;
        history_list_1{num_history_list_1,2}=tarlist_ofc;
        for k=1:num_t_ofc
            num2=1;
            num=num*num2;
        end
        cTOT=TOT/num_t_ofc;
        cVOC=VOC(c)/num_t_ofc;
        number_iter_2=1;
        len_2_min=-1000;                                 
        T2=com_parameter.KT2*num_t_ofc;
        number_iter_2max=ceil(log(0.001/T2)/log(com_parameter.K2));
        len_2=ones(1,1)*len_2_min;
        flag_2=zeros(1,1);
        flag_youjie=0;
        flag2=0;
        sumtc2=0;
        
        
        costmin=inf;
        i_cho=0;
        for k=1:num_shixu_c
            cost_sum=0;
            city_2=bh2list(tarlist_ofc,num_t_ofc,k);
            rvc_c=RVC;
            rvc_t=COE2RV(COET(:,city_2(1)));
            hc = cross(rvc_c(1:3), rvc_c(4:6));
            ht = cross(rvc_t(1:3), rvc_t(4:6));
            numda=COEC(1,c)/COET(1,city_2(1));
            thetaym=acos_2(dot(hc,ht)/(norm(hc)*norm(ht)));
            [theta1,theta2]=get2points(COEC(:,c),COET(:,city_2(1)));
            thetaph=theta2-theta1;
            if(thetaph<-pi)
                thetaph=thetaph+2*pi;
            elseif(thetaph>pi)
                thetaph=thetaph-2*pi;
            end
            cost=2*abs(sqrt(2-(1+thetaph/2/pi).^(-2/3))-1)+2*sin(thetaym/2)+abs(sqrt(2/(numda*numda+numda))-sqrt(1/numda))+abs(1-sqrt(2*numda/(numda+1)));
            cost_sum=cost_sum+cost;
            for j=2:num_t_ofc
                rvc_t_last=COE2RV(COET(:,city_2(j)));
                ht_last = cross(rvc_t_last(1:3), rvc_t_last(4:6));
                ht=ht_last;
                numda=COET(1,city_2(j-1))/COET(1,city_2(j));
                thetaym=acos_2(dot(ht,ht_last)/(norm(ht)*norm(ht_last)));
                [theta1,theta2]=get2points(COET(:,city_2(j-1)),COET(:,city_2(j)));
                thetaph=theta2-theta1;
                if(thetaph<-pi)
                    thetaph=thetaph+2*pi;
                elseif(thetaph>pi)
                    thetaph=thetaph-2*pi;
                end
                cost=2*abs(sqrt(2-(1+thetaph/2/pi).^(-2/3))-1)+2*sin(thetaym/2)+abs(sqrt(2/(numda*numda+numda))-sqrt(1/numda))+abs(1-sqrt(2*numda/(numda+1)));
                cost_sum=cost_sum+cost;
                
            end
            if(cost_sum<costmin)
                costmin=cost_sum;
                i_cho=k;
            end
        end
        
        
        city_2=bh2list(tarlist_ofc,num_t_ofc,i_cho);
        [len_2_min,flag_youjie,sumtc2,finalc]=func2(RTOM,m,Mu,c,num,city_2,num_t_ofc,RVC,cTOT,cVOC,COEC,COET,com_parameter);
        if(len_2_min<VOC(c))
            flag_youjie=1;
        else
            flag_youjie=0;
        end
        history_list_1{num_history_list_1,3}=len_2_min;   
        history_list_1{num_history_list_1,4}=sumtc2;      
        history_list_1{num_history_list_1,5}=flag_youjie;    
        history_list_1{num_history_list_1,6}=finalc;   
        
        if(flag_youjie==1)
            VOC_after(c)=VOC(c)-len_2_min;
            costfuel=costfuel+len_2_min;
            if(costfuel==0)
                a=1;
            end
            if(sumtc2>costtime)
                costtime=sumtc2;
            end
            final_pre2{c,1}=finalc;
        else
            flagi=0;
            break;
        end
    end
end
if(flagi)
    for u=1:m
        if(city(u)~=0)
            num_mi=num_mi+1;
        end
    end
    len=num_mi*2000-costfuel;
    sumt=costtime;
    sumfuel=costfuel;
    flag=1;
    final=final_pre2;
else
    len=-1;
    sumt=0;
    sumfuel=0;
    flag=0;
    final=final_pre2;
end
return

function [len,flagc,sumtc,finalc]=func2(RTOM,m,Mu,c,num,city_2,num_t_ofc,RVC,cTOT,cVOC,COEC,COET,com_parameter)
len=-1000;
flagc=0;
sumtc=0;
finalc=[];
city_3=ones(num_t_ofc,1);
if(isempty(city_2))
    return
end
[len_3,flag_3,sumtc_3,finalc_3]=func3(c,city_2,num_t_ofc,RTOM,cTOT,cVOC,COEC(:,c),COET);
if((len_3<len || len<-1)&& flag_3==1)
    flagc=1;
    finalc=finalc_3;
    len=len_3;
    sumtc=sumtc_3;
end
return

function [len,flagt,sumtc,finalc]=func3(c,city_2,num_t_ofc,RTOM,cTOT,cVOC,coec,COET)
RVC=COE2RV(coec);
RC=RVC(1:3);
VC=RVC(4:6);
t1=RTOM(3,c);                                 
pson=200;
psot=50;
psoc1=2;
psoc2=2;
psow1=0.95;
psow2=0.4;
psox=zeros(pson,num_t_ofc);
aa=zeros(1,100);
bb=zeros(num_t_ofc,100);
for k=1:1
    for u=1:num_t_ofc
        psox(:,u)=(cTOT*3600-t1)*ones(pson,1).*rand(pson,1);
    end
    vmax=0.75*(cTOT*3600-t1);
    v=rand(pson,num_t_ofc)*2*vmax-vmax;
    psop=psox;
    pbest=ones(pson,num_t_ofc);
    ppp=ones(pson,1);
    flag=ones(pson,1);
    ppp_best=inf;
    for i=1:pson
        [ppp(i),pbest(i,:),~,flag(i),~,~]=func13(c,RTOM,cTOT,cVOC,psox(i,:),RC,VC,coec,COET,num_t_ofc,city_2);
        if(ppp(i)<ppp_best&&flag(i)>0)
            ppp_best=ppp(i);
            psop_best=psox(i,:);
        end
    end
    psog=psop(1,:);
    gbest=inf;
    for i=1:pson
        if(ppp(i)<gbest && flag(i)==1)
            psog=psop(i,:);
            gbest=ppp(i);
        end
    end
    gb=ones(1,psot);
    flagyuejie=zeros(pson,psot);
    flagyuejie2=zeros(pson,psot);
    a=zeros(1,psot);
    cout_dead=0;
    i=0;
    while(i<psot || cout_dead<50)
        i=i+1;
        psow=psow1-(psow1-psow2)*i/psot;
        if(psot<i)
            psow=psow1;
        end
        for j=1:pson
            v(j,:)=psow*v(j,:)+psoc1*rand*(psop(j,:)-psox(j,:))+psoc2*rand*(psog-psox(j,:));
            psox(j,:)=psox(j,:)+v(j,:);
            flag_x=0;
            for ii=1:num_t_ofc
                if(v(j,ii)>vmax)
                    flagyuejie2(j,i)=1;
                    v(j,ii)=vmax;
                end
                if(v(j,ii)<-vmax)
                    flagyuejie2(j,i)=1;
                    v(j,ii)=-vmax;
                end
                if(psox(j,ii)<0  || (psox(j,ii)+t1)>cTOT*3600)
                    psox(j,ii)=(cTOT*3600-t1)*rand;
                    flagyuejie(j,i)=1;
                    flag_x=1;
                end
            end
            [ppp_2,pbest_2,~,flagt_2,~,~]=func13(c,RTOM,cTOT,cVOC,psox(j,:),RC,VC,coec,COET,num_t_ofc,city_2);
            if(flag(j)==0)
                flag(j)=flagt_2;
                psop(j,:)=psox(j,:);
                ppp(j)=ppp_2;
            elseif(flag(j)==1 && flagt_2==1)
                psop(j,:)=psox(j,:);
                ppp(j)=ppp_2;
            end
            if(flag(j)==1)
                if(ppp(j)<gbest)
                    psog=psop(j,:);
                    gbest=ppp(j);
                end
            end
        end
        if(gbest>1e10)
            break;
        end
        % a(i)=sum(flag);
        gb(i)=sum(gbest);
        if(i>1)
            cc=(gb(i-1)-gb(i))/gb(i-1);
            if(cc<0.01)
                cout_dead=cout_dead+1;
            else
                cout_dead=0;
            end
        end
    end
    aa(k)=gb(psot);
    bb(:,k)=psog;
    [len,len2,cost,flagt,sumtc,finalc]=func13(c,RTOM,cTOT,cVOC,psog,RC,VC,coec,COET,num_t_ofc,city_2);
end
return

function [len,len2,costk,flagt,sumtc,finalc]=func13(c,RTOM,cTOT,cVOC,x,RC,VC,coec,COET,num_t_ofc,city_2)

Mu=3.98600441500e14;    %miu,
flagt=1;
costk=0;
sumt=0;
revu=zeros(1,num_t_ofc);                                   
timeu1=zeros(1,num_t_ofc);                              
timeu2=zeros(1,num_t_ofc);                                  
timeu3=zeros(1,num_t_ofc);                               
deltavu=zeros(1,num_t_ofc);                               
deltav1=zeros(3,num_t_ofc);                               
deltav2=zeros(3,num_t_ofc);                              
dvgz=5;
for u=1:num_t_ofc 
    timeu1(u)=sumt;
    coet=COET(:,city_2(u));
    nT=sqrt(Mu/COET(1,city_2(u))^3);
    nTc=sqrt(8*Mu/(coec(1)+COET(1,city_2(u)))^3);
    nc=sqrt(Mu/coec(1)^3);
    tcc_cycle=2*pi/nc;
    tc_cycle=2*pi/nTc;
    coet(6)=coet(6)+nT*sumt;
    sumt=sumt+x(u);
    timeu2(u)=sumt;
    sumt=sumt+RTOM(3,c);
    timeu3(u)=sumt;
    rvc_after=COE2RV(coec);
    coet(6)=coet(6)+nT*x(u);
    rvt_after=COE2RV(coet);
    r1 = norm(rvc_after(1:3));
    r2 = norm(rvt_after(1:3));
    c12 = cross(rvc_after(1:3)/r1, rvt_after(1:3)/r2);
    dtheta = acos_2(dot(rvc_after(1:3)/r1,rvt_after(1:3)/r2));
    
    if c12(3) <= 0
        dtheta = 2*pi - dtheta;
    end
    phi_zc=x(u)/tcc_cycle-floor(x(u)/tcc_cycle);
    rev=round(x(u)/tc_cycle);
    if(phi_zc< pi)
        if(dtheta > (pi+phi_zc))
            rev=rev-1;
            if(rev<0)
                rev=0;
            end
        end
    else
        if(dtheta > (phi_zc-pi))
            rev=rev-1;
        end
    end
    VC=rvc_after(4:6);
    [V1,V2,flag]=lambertu(VC,rvc_after(1:3),rvt_after(1:3),x(u),'pro',rev);
    VT= rvt_after(4:6);
    num_v=1;
    if(flag==2)
        [mindv,num_v]=min([norm(V1(:,1)-VC)+norm(V2(:,1)-VT),norm(V1(:,2)-VC)+norm(V2(:,2)-VT)]);
    elseif(flag==1)
        mindv=norm(V1-VC)+norm(V2-VT);
    else
        flagt=0;
        revu(u)=rev+1;
        break;
    end
    coet(6)=coet(6)+nT*RTOM(3,c);
    coec=coet;
    mindv=mindv+dvgz;
    revu(u)=rev+1;
    deltavu(u)=mindv;
    deltav1(:,u)=V1(:,num_v)-VC;
    deltav2(:,u)=VT-V2(:,num_v);
    costk=costk+mindv;
    len2(u)=mindv;
end
len=costk;
sumtc=sumt;
if(flagt>0)
    finalc=zeros(16,1);
    for u=1:num_t_ofc 
        finalc(1,u)=city_2(u);
        finalc(2,u)=revu(u);
        finalc(3,u)=deltavu(u);
        finalc(4,u)=timeu1(u);
        finalc(5,u)=timeu2(u);
        finalc(6,u)=timeu3(u);
        finalc(7,u)=deltav1(1,u);
        finalc(8,u)=deltav1(2,u);
        finalc(9,u)=deltav1(3,u);
        finalc(10,u)=deltav2(1,u);
        finalc(11,u)=deltav2(2,u);
        finalc(12,u)=deltav2(3,u);
    end
else
    for u=1:num_t_ofc 
        len2(u)=0;
    end
    flagt=0;
    finalc=[];
end 
return

function bb=acos_2(aaa)
if(aaa>=1)
    aaa=1;
elseif(aaa<=-1)
    aaa=-1;
end
bb=acos(aaa);
return


function kdv=cal_k(n,m,TOT)
Mu=3.98600441500e14;    %miu,开普勒常数*4*pi*pi对应的值，单位为m^3/s^2
timeb=TOT*3600/m*n;
coec=[42164169;0;0;0;0;0];
coet=[42164169;0;0;0;0;deg2rad(30)];
coet2=[42164169;0;0;0;0;deg2rad(-30)];
nT=sqrt(Mu/coet(1)^3);
TTT=2*pi/nT;
if(timeb>TTT)
    timeb=TTT*floor(timeb/TTT);
end
RVC=COE2RV(coec);
RC=RVC(1:3);
VC=RVC(4:6);
rev=floor(timeb/TTT)-1;
if(rev<0)
    rev=0;
end
timeb=timeb-TTT/12;
coet(6)=coet(6)+nT*timeb;
RVT=COE2RV(coet);
RT=RVT(1:3);
VT=RVT(4:6);
    [V1,V2,flag]=lambertu(VC,RC,RT,timeb,'pro',rev);
    [V12,V22,flag2]=lambertu(VC,RC,RT,timeb,'pro',rev+1);
    [V13,V23,flag3]=lambertu(VC,RC,RT,timeb,'pro',rev+2);
    deltav=[];
    if (flag~=0)
        numv=size(V1,2);
        if(numv==2)
            DeltaV11=V1(:,1)-VC;
            DeltaV21=VT-V2(:,1);
            DeltaV12=V1(:,2)-VC;
            DeltaV22=VT-V2(:,2);
            v1s=[norm(DeltaV11)+norm(DeltaV21),norm(DeltaV12)+norm(DeltaV22)];
            [deltav1,~]=min(v1s);
        else
            DeltaV1=V1-VC;
            DeltaV2=VT-V2;
            deltav1=norm(DeltaV1)+norm(DeltaV2);
        end
        deltav=[deltav,deltav1];
    end
    if (flag2~=0)
        numv=size(V12,2);
        if(numv==2)
            DeltaV11=V12(:,1)-VC;
            DeltaV21=VT-V22(:,1);
            DeltaV12=V12(:,2)-VC;
            DeltaV22=VT-V22(:,2);
            v1s=[norm(DeltaV11)+norm(DeltaV21),norm(DeltaV12)+norm(DeltaV22)];
            [deltav2,~]=min(v1s);
        else
            DeltaV1=V12-VC;
            DeltaV2=VT-V22;
            deltav2=norm(DeltaV1)+norm(DeltaV2);
        end
        deltav=[deltav,deltav2];
    end
    if (flag3~=0)
        numv=size(V13,2);
        if(numv==2)
            DeltaV11=V13(:,1)-VC;
            DeltaV21=VT-V23(:,1);
            DeltaV12=V13(:,2)-VC;
            DeltaV22=VT-V23(:,2);
            v1s=[norm(DeltaV11)+norm(DeltaV21),norm(DeltaV12)+norm(DeltaV22)];
            [deltav3,~]=min(v1s);
        else
            DeltaV1=V13-VC;
            DeltaV2=VT-V23;
            deltav3=norm(DeltaV1)+norm(DeltaV2);
        end
        deltav=[deltav,deltav3];
    end
    
    deltav1111=min(deltav);%+t/400;
    
timeb=timeb+TTT/6;
coet2(6)=coet2(6)+nT*timeb;
RVT2=COE2RV(coet2);
RT2=RVT2(1:3);
VT2=RVT2(4:6);
    [V1,V2,flag]=lambertu(VC,RC,RT2,timeb,'pro',rev);
    [V12,V22,flag2]=lambertu(VC,RC,RT2,timeb,'pro',rev+1);
    [V13,V23,flag3]=lambertu(VC,RC,RT2,timeb,'pro',rev+2);
    deltav=[];
    if (flag~=0)
        numv=size(V1,2);
        if(numv==2)
            DeltaV11=V1(:,1)-VC;
            DeltaV21=VT2-V2(:,1);
            DeltaV12=V1(:,2)-VC;
            DeltaV22=VT2-V2(:,2);
            v1s=[norm(DeltaV11)+norm(DeltaV21),norm(DeltaV12)+norm(DeltaV22)];
            [deltav1,~]=min(v1s);
        else
            DeltaV1=V1-VC;
            DeltaV2=VT2-V2;
            deltav1=norm(DeltaV1)+norm(DeltaV2);
        end
        deltav=[deltav,deltav1];
    end
    if (flag2~=0)
        numv=size(V12,2);
        if(numv==2)
            DeltaV11=V12(:,1)-VC;
            DeltaV21=VT2-V22(:,1);
            DeltaV12=V12(:,2)-VC;
            DeltaV22=VT2-V22(:,2);
            v1s=[norm(DeltaV11)+norm(DeltaV21),norm(DeltaV12)+norm(DeltaV22)];
            [deltav2,~]=min(v1s);
        else
            DeltaV1=V12-VC;
            DeltaV2=VT2-V22;
            deltav2=norm(DeltaV1)+norm(DeltaV2);
        end
        deltav=[deltav,deltav2];
    end
    if (flag3~=0)
        numv=size(V13,2);
        if(numv==2)
            DeltaV11=V13(:,1)-VC;
            DeltaV21=VT2-V23(:,1);
            DeltaV12=V13(:,2)-VC;
            DeltaV22=VT2-V23(:,2);
            v1s=[norm(DeltaV11)+norm(DeltaV21),norm(DeltaV12)+norm(DeltaV22)];
            [deltav3,~]=min(v1s);
        else
            DeltaV1=V13-VC;
            DeltaV2=VT2-V23;
            deltav3=norm(DeltaV1)+norm(DeltaV2);
        end
        deltav=[deltav,deltav3];
    end
    
    deltav2222=min(deltav);%+t/400;
    kdv=(deltav1111+deltav2222)/60;
return

function [list_serving_codeold,list_serving_code,one2one,num_serving_T,num_fenpei]=onetoone(n,m,RTOM,VOC,TOT,COEC,COET)
list_serving_codeold=cell(m,1);                             
num_t_ofc=1;
num_fenpei=1;
list_serving_code=zeros(m,n);
for i=1:m
    p=0;
    for c=1:n
        cTOT=TOT;
        cVOC=VOC(c);
        city_2=i;
        coec=COEC(:,c);
       [len,flagt,sumtc,finalc]=func3(c,city_2,num_t_ofc,RTOM,cTOT,cVOC,coec,COET);
       if(flagt)
           p=p+1;
           list_serving_code(i,c)=1;
           list_serving_codeold{i,1}(1,p)=c;
           one2one{i,c}=finalc;
       else
           one2one{i,c}=-1;
       end
    end
    num_serving_T(i)=p;
    num_fenpei=num_fenpei*p;
end

return

function [theta1,theta2] = get2points( coe1,coe2 ) %(0-2pi]
rvc1=COE2RV(coe1);
rvc2=COE2RV(coe2);
h1=cross(rvc1(1:3),rvc1(4:6));
h1=h1/norm(h1);
h2=cross(rvc2(1:3),rvc2(4:6));
h2=h2/norm(h2);
line=cross(h1,h2);
if (all(line==0))
    theta1=2*pi;
    theta2=2*pi;
else
    c11 = cross(line,rvc1(1:3));
    c12 = cross(line,rvc2(1:3));
    theta1=acos(dot(line,rvc1(1:3))/(norm(line)*norm(rvc1(1:3))));
    theta2=acos(dot(line,rvc2(1:3))/(norm(line)*norm(rvc2(1:3))));
    if c11(3)>= 0
        theta1 = 2*pi - theta1;
    end
    if c12(3)>= 0
        theta2 = 2*pi - theta2;
    end
end
return


function [V1, V2, Flag] = lambertu(V1_PRE,R1, R2, t, type, rev, mu)
Flag=1;
flag=0;
FLAG=0;
if nargin < 7
    mu = 3.98600436e14; 
end

if nargin < 6
  rev = 0;  
end

if nargin < 5
   type = 'pro';
end
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1,R2)/r1/r2);

if ~(strcmp(type, 'pro') || strcmp(type,'retro'))
    type = 'pro';
    fprintf('\n ** Prograde trajectory assumed.\n')
end

if strcmp(type, 'pro')
    if c12(3) <= 0
        theta = 2*pi - theta;
    end
else 
    if c12(3) >= 0
        theta = 2*pi - theta;
    end
end


A = sin(theta) * sqrt(r1 * r2 / (1 - cos(theta)));
if(theta<0.001 && theta>0)
    if(norm((r2-r1)/r1)<0.002)
        if(rev<1)
            rev=1;
        end
        t_tx=t/(rev);
        a_tx=(0.5*t_tx*sqrt(mu)/pi).^(2/3);
        T_zq=2*pi*sqrt(r1.^3/mu);
        if t_tx>T_zq
            e_tx=(a_tx-r1)/a_tx;
            if((1+e_tx)<0)
                V2 =0;
                V1=V2;
                Flag=0;
                return
            else
                h=sqrt(r1*mu*(1+e_tx));
            end
        else
            e_tx=(r1-a_tx)/a_tx;
            if((1-e_tx)<0)
                V2 =0;
                V1=V2;
                Flag=0;
                return
            else
                h=sqrt(r1*mu*(1-e_tx));
            end
        end
        v_1=h/r1;
        rate_v=v_1/norm(V1_PRE);
        V2 = V1_PRE*rate_v;
        V1=V2;
        Flag=1;
        return
    else
        V2 =0;
        V1=V2;
        Flag=0;
        return
    end
end
if(theta-2*pi<=0 && theta-2*pi>-0.001)
    if(norm((r2-r1)/r1)<0.002)
        t_tx=t/(rev+1);
        a_tx=(0.5*t_tx*sqrt(mu)/pi).^(2/3);
        T_zq=2*pi*sqrt(r1.^3/mu);
        if t_tx>T_zq
            e_tx=(a_tx-r1)/a_tx;
            if((1+e_tx)<0)
                V2 = 0;
                V1=V2;
                Flag=0;
                return
            else
                h=sqrt(r1*mu*(1+e_tx));
            end
        else
            e_tx=(r1-a_tx)/a_tx;
            if((1-e_tx)<0)
                V2 = 0;
                V1=V2;
                Flag=0;
                return
            else
                h=sqrt(r1*mu*(1-e_tx));
            end
        end
        v_1=h/r1;
        rate_v=v_1/norm(V1_PRE);
        V2 = V1_PRE*rate_v;
        V1=V2;
        Flag=1;
        return
    else
        V2 =0;
        V1=V2;
        Flag=0;
        return
    end
end
if(-1e-10<A && A<1e-10)
     V2 =0;
     V1=V2;
     Flag=0;
     return
end
if rev == 0
    z = [(2*pi).^2-0.001  35 30 25 20 15 10 5 1 0.001 0.000001];
    [Fz, yz] = F(z,r1,r2,A);
    ia = find(yz>0, 1, 'last');
    ib = find(yz<0, 1, 'first'); % z>=zmin, 
    if ~isempty(ib)
        % y(z) = 0
        tol = 1e-10;
        [FLAG,zmin] = bisection(@(x) y(x,r1,r2,A),z(ib), z(ia), tol );
        za = zmin;
    else
        ind = find(Fz<sqrt(mu)*t, 1, 'first');
        if isempty(ind)
            %printf('\n varable z is too small.');
            za = z(end);
        else
            za = z(ind);
        end
    end 
    zb = (2*pi)^2 - 0.001;
    tol = 1e-10;
   
    [flag,z] = bisection(@(x) (F(x,r1,r2,A) - sqrt(mu) * t), za, zb, tol );   
    if (flag==1)
        V1=0;
        V2=0;
        Flag=0;
        return
    end
    [Fz,yz] = F(z,r1,r2,A); 
    f = 1 - yz/r1;   
    g = A*sqrt(yz/mu);   
    gdot = 1 - yz/r2;
    V1 = 1/g*(R2 - f*R1);    
    V2 = 1/g*(gdot*R2 - R1);
    
    V1 = V1(:);
    V2 = V2(:); 
else
   
   ia =[];
   ib =[];
   z1 = (2*pi*(rev)).^2;
   z2 = (2*pi*(rev+1)).^2;
   while isempty(ia) || isempty(ib)
       z = linspace(z1,z2,101);
       dFz = dFdz(z,r1,r2,A);
       ia = find(dFz<0, 1, 'last');
       ib = find(dFz>0, 1, 'first');
       if ~isempty(ia)
           z1 = z(ia);
       end
       if ~isempty(ib)
           z2 = z(ib);
       end
   end
   za = z(ia);
   zb = z(ib);
   
   tol = 1e-10;
   [FLAG,zmin] = bisection(@(x) dFdz(x,r1,r2,A), za, zb, tol );
    if (FLAG==1)
        V1=0;
        V2=0;
        Flag=0;
        return
    end
   [Fzmin, yz]= F(zmin,r1,r2,A);
   if ((Fzmin-sqrt(mu) * t)/Fzmin<0.001 && (Fzmin-sqrt(mu) * t)/Fzmin>-0.001) 
       %...Equation 5.46a:
       f = 1 - yz/r1;
       %...Equation 5.46b:
       g = A*sqrt(yz/mu);
       %...Equation 5.46d:
       gdot = 1 - yz/r2;
       %...Equation 5.28:
       V1 = 1/g*(R2 - f*R1);
       %...Equation 5.29:
       V2 = 1/g*(gdot*R2 - R1);
       V1 = V1(:);
       V2 = V2(:);       
   elseif Fzmin < sqrt(mu) * t  
       % left part
       zb = zmin;
       za = (2*pi*(rev)).^2 + 0.001;
       tol = 1e-10;
       [FLAG_1,z] = bisection(@(x) (F(x,r1,r2,A) - sqrt(mu) * t), za, zb, tol );
       if (FLAG_1~=1)
           [Fz,yz] = F(z,r1,r2,A);
           %...Equation 5.46a:
           f = 1 - yz/r1;
           %...Equation 5.46b:
           g = A*sqrt(yz/mu);
           %...Equation 5.46d:
           gdot = 1 - yz/r2;
           %...Equation 5.28:
           V1L = 1/g*(R2 - f*R1);
           %...Equation 5.29:
           V2L = 1/g*(gdot*R2 - R1);
       end
       % right part
       za = zmin;
       zb = (2*pi*(rev+1)).^2 - 0.001;
       tol = 1e-10;
       [FLAG_2,z] = bisection(@(x) (F(x,r1,r2,A) - sqrt(mu) * t), za, zb, tol );
       if (FLAG_2~=1)
           [Fz,yz] = F(z,r1,r2,A);
           %...Equation 5.46a:
           f = 1 - yz/r1;
           %...Equation 5.46b:
           g = A*sqrt(yz/mu);
           %...Equation 5.46d:
           gdot = 1 - yz/r2;
           %...Equation 5.28:
           V1R = 1/g*(R2 - f*R1);
           %...Equation 5.29:
           V2R = 1/g*(gdot*R2 - R1);
       end
       if (FLAG_1&&FLAG_2)
           V1=0;
           V2=0;
           Flag=0;
       elseif(FLAG_1)
           V1=V1R(:);
           V2=V2R(:);
           Flag=1;
       elseif(FLAG_2)
           V1=V1L(:);
           V2=V2L(:);
           Flag=1;
       else
           V1 = [V1L(:),V1R(:)];
           V2 = [V2L(:),V2R(:)];
           Flag=2;
       end
       return
   else 
       V1=0;
       V2=0;
       Flag=0;      
   end
   
end

return



% Subfunctions used in the main body:
function [yz,cz,sz] = y(z,r1,r2,A)
cz = stumpC(z);
sz = stumpS(z);
yz = r1 + r2 + A*(z.*sz - 1) ./ sqrt(cz);
return

%...Equation 5.40:
function [Fz,yz] = F(z,r1,r2,A)
[yz,cz,sz] = y(z,r1,r2,A);
Fz = (yz./cz).^1.5 .* sz + A * sqrt(yz);
return

%...Equation 5.43:
function [dFz,yz] = dFdz(z,r1,r2,A)
[yz,cz,sz] = y(z,r1,r2,A);

dFz = (yz./cz).^1.5 .* (1/2 * (cz - 3/2* sz./cz) ./ z ...
    + 3*sz.^2/4./cz) ...
    + A/8*(3*sz./cz.*sqrt(yz) ...
    + A*sqrt(cz./yz));

yz0 = yz(z==0);
dFz(z==0) = sqrt(2)/40 * yz0.^1.5 + A/8 * (sqrt(yz0) ...
    + A * sqrt(1/2*ones(size(yz0))./yz0));
return

function [flag,x0] = bisection(func, a, b, tol )
% Bisection method for finding function root.

fa = func(a);
fb = func(b);


if fa * fb > 0
    x0=0;
    flag=1;
    return
end

maxItr = ceil((log(b-a)-log(tol)) / log(2.0));

for i=1:maxItr
    c = (a + b) / 2.0;
    fc = func(c);
    
    if abs(fc) < eps
        x0 = c;
        flag=0;
        return
    end
    
    if fb * fc > 0
        b = c;
        fb = fc;
    else
        a = c;
        fa = fc;
    end
    
    if (b-a) < tol
        x0 = (a+b)/2;
        flag=0;
        return
    end
end
flag=0;
x0 = c;
return

function X = COE2RV(coe)
if nargin<1
    coe=[42166e3,0.001,rad2deg(7),0,0,pi]';
end

aT=coe(1);
eT=coe(2);
iT=coe(3);
raanT=coe(4);
omegaT=coe(5);
fT=M2f(coe(6),coe(2));

Mu = 398600.436e+9;     

MZ1=[cos(-raanT),sin(-raanT),0.0;
       -sin(-raanT),cos(-raanT),0.0;
       0.0,0.0,1.0];
MX=[1.0,0.0,0.0;
       0.0,cos(-iT),sin(-iT);
       0.0,-sin(-iT),cos(-iT)];
MZ2=[cos(-(omegaT + fT)),sin(-(omegaT + fT)),0.0;
       -sin(-(omegaT + fT)),cos(-(omegaT + fT)),0.0;
       0.0,0.0,1.0];
TR=MZ1*MX*MZ2;
r_t = aT * (1.0 - eT^2) / (1.0 + eT * cos(fT));
vt_x = sqrt(Mu / (aT * (1.0 - eT^2))) * eT * sin(fT);
vt_y = sqrt(Mu / (aT * (1.0 - eT^2))) * (1.0 + eT * cos(fT));
X(1:3,1) = TR * [r_t;0.0;0.0];
X(4:6,1) = TR * [vt_x;vt_y;0.0];
return

function f = M2f(M,e)
if nargin<2
    e=0.01;
    M=2*pi;
end
M=atan2(sin(M),cos(M));
if M<0
    M=M+2*pi;
end
if M<=pi
    E1=0;
    E2=pi;
    for k=1:100000
        E=(E1+E2)/2;
        if abs(E2-E1)<1e-6
            break;
        else
            M3=E-e*sin(E);
            if M3>M
                E2=E;
            else
                E1=E;
            end
        end
    end
else
    E1=pi;
    E2=2*pi;
    for k=1:100000
        E=(E1+E2)/2;
        if abs(E2-E1)<1e-6
            break;
        else
            M3=E-e*sin(E);
            if M3>M
                E2=E;
            else
                E1=E;
            end
        end
    end
end
cosf=(cos(E)-e)/(1-e*cos(E));
sinf=sqrt(1-e^2)*sin(E)/(1-e*cos(E));
f=atan2(sinf,cosf);
if f<0
    f=f+pi*2;
end
return

function c = stumpC(z)



c = zeros(size(z));

ze = z(z > 0);
zh = z(z < 0); 

% z > 0, ellipse 
c(z > 0) = (1 - cos(sqrt(ze))) ./ ze;

% z == 0, parabola
c(z == 0) = 1/2;

% z < 0, hyperbola
c(z < 0) = (cosh(sqrt(-zh)) - 1) ./ (-zh);

return

function s = stumpS(z)


s = zeros(size(z));

ze = z(z > 0);
zh = z(z < 0); 

% z > 0, ellipse 
s(z > 0) = (sqrt(ze) - sin(sqrt(ze))) ./ (sqrt(ze)).^3;

% z == 0, parabola
s(z == 0) = 1/6;

% z < 0, hyperbola
s(z < 0) = (sinh(sqrt(-zh)) - sqrt(-zh)) ./ (sqrt(-zh)).^3;

return

     
