function data_mean=load_pu(test_class)
load paviaU_gt
load paviaU


data_in=paviaU;
gt=paviaU_gt;
for i=1:size(data_in,3)
    a=data_in(:,:,i);
    data_vec(:,i)=a(:);
end
data=[];
class=[];
for i=1:numel(unique(gt))-1
    locs=find(gt==i);
    data=[data;data_vec(locs,:)];
    class=[class;ones(numel(locs),1)*i];
end
total_data=data;
total_class=class;

for i=1:numel(test_class)
    data=total_data(total_class==test_class(i),:);
    data_mean(i,:)=mean(data);
end