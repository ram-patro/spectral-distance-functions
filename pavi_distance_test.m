clear
close all
clc

prob=[0.1 0.02 0.05 0.75 0.08];

%Pavia University
test_class=[2 3 6 7 8];
data_mean=load_pu(test_class);
plot(data_mean','LineWidth',2);xlim([1 103]);legend('Meadows','Gravel','Bare Soil','Bitumen','Self-Blocking Bricks');grid on;
xlabel('Bands');ylabel('Reflectance');

classes=numel(test_class);

D=sum(data_mean.*prob');

method={'CBD','ED','NED','DPD','SAM','SCA','CCSM','SSV','SID','NS3','SGA','SID_SAM_tan','SID_SAM_sin','SID_SCA_tan','SID_SCA_sin'}';

for k=1:numel(method)
    m{1}=method{k};
    DIST.(method{k}).d=zeros(classes,classes);
    for i=1:classes-1
        a=data_mean(i,:);
        for j=i+1:classes
            b=data_mean(j,:);
            DIST.(method{k}).d(i,j)=spectral_distance_all(a,b,m);
        end
    end
end


for k=1:numel(method)
    m{1}=method{k};
    SDPW.(method{k}).d=zeros(classes,classes);
    for i=1:size(data_mean,1)-1
        a=data_mean(i,:);
        for j=i+1:size(data_mean,1)
            b=data_mean(j,:);
            SDPW.(method{k}).d(i,j)=calc_SDPW(a,b,D,m);
        end
    end
end

for k=1:numel(method)
    m{1}=method{k};
    SDPR.(method{k}).d=zeros(1,classes);
    SDEN.(method{k}).d=0;
    for i=1:classes
        a=data_mean(i,:);
        sum_all=zeros(1,classes);
        for j=1:classes
            b=data_mean(j,:);
            sum_all(j)=spectral_distance_all(b,D,m);
        end
        SDPR.(method{k}).d(i)=spectral_distance_all(a,D,m)/sum(sum_all);
    end
    SDEN.(method{k}).d=-sum((SDPR.(method{k}).d).*log(SDPR.(method{k}).d));
end

save DIST DIST
save SDPW SDPW
save SDPR SDPR
save SDEN SDEN
save method method
save prob prob
