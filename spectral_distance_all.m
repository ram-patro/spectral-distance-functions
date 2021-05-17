function d=spectral_distance_all(a,b,type)
global time
d_all=zeros(1,2);
tic;
for i=1:numel(type)
    if(strcmp(type{i},'CBD'))%City block distance
        d=hyper_CBD(a,b);
    elseif(strcmp(type{i},'ED'))%Euclidian    
        d=hyper_ED(a,b);
    elseif(strcmp(type{i},'NED'))%Normalized Euclidian    
        d=hyper_NED(a,b);    
    elseif(strcmp(type{i},'DPD'))%Dot Product distance    
        d=hyper_DPD(a,b);   
    elseif(strcmp(type{i},'SAM'))%Spectral angle Mapper
       d=hyper_SAM(a,b);
    elseif(strcmp(type{i},'SCA'))%Spectral corelation Mapper
       d=hyper_SCA(a,b);
    elseif(strcmp(type{i},'CCSM'))%Cross correlogram spectral matcher
       d=hyper_CCSM(a,b);
    elseif(strcmp(type{i},'SSV'))%Spectral Similarity Value
        d=hyper_SSV(a,b);
    elseif(strcmp(type{i},'SID'))%Spectral Information Divergence
        d=hyper_SID(a,b);
    elseif(strcmp(type{i},'NS3'))%Normalized spectral similarity score
        d=hyper_NS3(a,b);    
    elseif(strcmp(type{i},'SGA'))%Spectral Graient Angle
        d=hyper_SGA(a,b);    
    elseif(strcmp(type{i},'SID_SAM_tan')) %SID-SAM-tan  
        d=hyper_SIDSAM_tan(a,b);
    elseif(strcmp(type{i},'SID_SAM_sin')) %SID-SAM-sin  
        d=hyper_SIDSAM_sin(a,b);
    elseif(strcmp(type{i},'SID_SCA_tan')) %SID-SCA-tan  
        d=hyper_SIDSCA_tan(a,b);
    elseif(strcmp(type{i},'SID_SCA_sin')) %SID-SCA-sin  
        d=hyper_SIDSCA_sin(a,b);
    end
    time=toc;
    d_all(i)=d;
end
end

function d=hyper_CBD(a,b)
    d=sum(abs(a-b));
end

function d=hyper_ED(a,b)
    d=sqrt(sum((a-b).^2));
end

function d=hyper_NED(a,b)
    diff=abs(a-b);
    d_min=min(diff);d_max=max(diff);
    dc=sqrt(mean((diff).^2));
    d=(dc-d_min)/(d_max-d_min);
end

function d=hyper_DPD(a,b)
    d=(sum(a.*b).^2)/((sum(a)^2)*(sum(b)^2));
end

function d=hyper_SAM(a,b)
    d=acos(dot(a, b)/ (norm(b) * norm(a)));
end

function d=hyper_SCA(a,b)
    r=dot(a-mean(a), b-mean(b))/ (norm(b-mean(b)) * norm(a-mean(b)));
    d=acos((r+1)/2);%radians
end

function d=hyper_CCSM(a,b)
cnt=1;
    for m=-10:10
        if(m<0)
            aa=a(1:end+m+1);
            bb=b(abs(m):end);
        elseif(m==0)
            aa=a;bb=b;
        else
            aa=a(abs(m):end);
            bb=b(1:end-m+1);
        end
        n=numel(aa);
        rm(cnt)=(n*sum(aa.*bb) -(sum(aa)*sum(bb)))/sqrt((n*sum(aa.^2)-sum(aa)^2)*(n*sum(bb.^2)-sum(bb)^2));
        cnt=cnt+1;
    end
    %plot(-10:10,rm)
    d=1-max(rm);
end

function d=hyper_SSV(a,b)
    n=numel(a);
    r2=((1/(n-1))*sum((a-mean(a)).*(b-mean(b))))/(std(a)*std(b));
    d=sqrt(hyper_ED(a,b)+(1-r2));
end

function d=hyper_SID(a,b)
    a=a/sum(a);b=b/sum(b);
    d=abs(sum(a.*log(a./b)) + sum(b.*log(b./a)));
end

function d=hyper_NS3(a,b)
    d=sqrt(hyper_NED(a,b)^2+(1-cos(hyper_SAM(a,b))^2));
end

function d=hyper_SGA(a,b)
    d=hyper_SAM([a(1:end-1)-a(2:end)],[b(1:end-1)-b(2:end)]);
end

function d=hyper_SIDSAM_sin(a,b)
    d=hyper_SID(a,b)*sin(hyper_SAM(a,b));    
end

function d=hyper_SIDSAM_tan(a,b)
    d=hyper_SID(a,b)*tan(hyper_SAM(a,b));
end

function d=hyper_SIDSCA_sin(a,b)
    d=hyper_SID(a,b)*sin(hyper_SCA(a,b));    
end

function d=hyper_SIDSCA_tan(a,b)
    d=hyper_SID(a,b)*tan(hyper_SCA(a,b));
end

