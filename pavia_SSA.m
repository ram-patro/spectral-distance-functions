clear
close all
clc

load paviaU
[R,C,B]=size(paviaU);
data_Vec=reshape(paviaU,[R*C B]);
M=24;
L=5;
no=17000;%arbitrary pixel value [1 - 207400]
[datassa(no,:),~,~]=ssa_recon(data_Vec(no,:),M,1:L);
plot(data_Vec(no,:),'LineWidth',2);hold on;
plot(datassa(no,:),'LineWidth',2);hold off;
xlim([1 103]);legend('Raw spectrum','SSA Reconstructed spectrum');
xlabel('Bands');ylabel('Reflectance');
