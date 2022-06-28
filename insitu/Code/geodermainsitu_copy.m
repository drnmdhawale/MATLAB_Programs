%%
%% A program to load the in-lab measurements on 20 locations(soil samples)
clear all
clc
exsitu1  = xlsread('VIS_NIR_In_Situ.xlsx', 'ALL_RW');exsitu1 =exsitu1';
exsitu2  = xlsread('VIS_NIR_In_Situ.xlsx', 'ALL_F');exsitu2 =exsitu2';
exsitu3  = xlsread('VIS_NIR_In_Situ.xlsx', 'ALL_FD');exsitu3 =exsitu3';
exsitu4  = xlsread('VIS_NIR_In_Situ.xlsx', 'ALL_SD');exsitu4 =exsitu4';

%% Combine factors
vis_nir_spectra=[exsitu1(6:end,2:end); exsitu2(6:end,2:end); exsitu3(6:end,2:end); ...
    exsitu4(6:end,2:end); exsitu1(2:5,2:end)];

vis_nir_factors=[exsitu1(6:end,1); exsitu2(6:end,1); ...
    exsitu3(6:end,1); exsitu4(6:end,1);];

wavelengths=vis_nir_factors(1:363,:);
%% calculate average spectrum
j=1;
for i=1:1:19,
    avg1(:,j)=nanmean([vis_nir_spectra(:,i) vis_nir_spectra(:,i+19) vis_nir_spectra(:,i+38)],2);
    var1(:,j)=nanvar([vis_nir_spectra(:,i) vis_nir_spectra(:,i+19) vis_nir_spectra(:,i+38)]')';
    j=j+1;
end

%%  Calculate stats
den= sqrt(mean(var1,2));
num=std(avg1')';
RSE=num./den;
VAR_DATA2=var(vis_nir_spectra')';
MST=var(avg1')'*3;
MSE=mean(var1,2);
TOTAL=VAR_DATA2*56;
SST=MST*19;
SSE=MSE*38;
F_MEAN=MST./MSE;
F_IND=F_MEAN/3;
RMSE=power(MSE, 1/2);
%RSE2=power(F_MEAN/3,1/2);

figure(1);
plot(vis_nir_factors(1:363,:), RSE(1:363,:), 'or'); hold on
plot(vis_nir_factors(364:726,:), RSE(364:726,:), 'xb'); hold on
plot(vis_nir_factors(727:1089,:), RSE(727:1089,:), '*k'); hold on
plot(vis_nir_factors(1090:1452,:), RSE(1090:1452,:), '+m'); grid off
xlabel('Wavelengths, nm', 'FontSize',16)
ylabel('RSE ', 'FontSize',16)
set(gca,'XTick',400:200:2200); set(gca,'YTick',0.0:1:12);
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','NorthEast')
xlim([400 2200])
ylim([0.0 12.0])
 
%% calculate correlations with Raw Spectrum

ll=1; hl=363;
Y=vis_nir_spectra(1453:1456,1:19)';
X=vis_nir_spectra(ll:hl,1:19)';


for i=1:1:363
        RHOt11(i,1) = corr(X(:,i),Y(:,1));
        RHOt12(i,1) = corr(X(:,i),Y(:,2));
        RHOt13(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,20:38)';

for i=1:1:363
        RHOt21(i,1) = corr(X(:,i),Y(:,1));
        RHOt22(i,1) = corr(X(:,i),Y(:,2));
        RHOt23(i,1) = corr(X(:,i),Y(:,3));
end
  
X=vis_nir_spectra(ll:hl,39:57)';

for i=1:1:363
        RHOt31(i,1) = corr(X(:,i),Y(:,1));
        RHOt32(i,1) = corr(X(:,i),Y(:,2));
        RHOt33(i,1) = corr(X(:,i),Y(:,3));
end

X=avg1(ll:hl,:)';
for i=1:1:363
        RHOM1(i,1) = corr(X(:,i),Y(:,1));
        RHOM2(i,1) = corr(X(:,i),Y(:,2));
        RHOM3(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,:)';
Y=vis_nir_spectra(1453:1456,:)';

for i=1:1:363
        RHOall1(i,1) = corr(X(:,i),Y(:,1));
        RHOall2(i,1) = corr(X(:,i),Y(:,2));
        RHOall3(i,1) = corr(X(:,i),Y(:,3));
end

corr_raw= [RHOall1, RHOall2 RHOall3, RHOM1, RHOM2, RHOM3, RHOt11, RHOt12, RHOt13, ...
    RHOt21, RHOt22, RHOt23, RHOt31, RHOt32, RHOt33];

%% calculate correlations with Filtered Spectrum
ll=364; hl=726;
Y=vis_nir_spectra(1453:1456,1:19)';
X=vis_nir_spectra(ll:hl,1:19)';


for i=1:1:363
        RHOt11(i,1) = corr(X(:,i),Y(:,1));
        RHOt12(i,1) = corr(X(:,i),Y(:,2));
        RHOt13(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,20:38)';

for i=1:1:363
        RHOt21(i,1) = corr(X(:,i),Y(:,1));
        RHOt22(i,1) = corr(X(:,i),Y(:,2));
        RHOt23(i,1) = corr(X(:,i),Y(:,3));
end
  
X=vis_nir_spectra(ll:hl,39:57)';

for i=1:1:363
        RHOt31(i,1) = corr(X(:,i),Y(:,1));
        RHOt32(i,1) = corr(X(:,i),Y(:,2));
        RHOt33(i,1) = corr(X(:,i),Y(:,3));
end

X=avg1(ll:hl,:)';
for i=1:1:363
        RHOM1(i,1) = corr(X(:,i),Y(:,1));
        RHOM2(i,1) = corr(X(:,i),Y(:,2));
        RHOM3(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,:)';
Y=vis_nir_spectra(1453:1456,:)';

for i=1:1:363
        RHOall1(i,1) = corr(X(:,i),Y(:,1));
        RHOall2(i,1) = corr(X(:,i),Y(:,2));
        RHOall3(i,1) = corr(X(:,i),Y(:,3));
end

corr_filt= [RHOall1, RHOall2 RHOall3, RHOM1, RHOM2, RHOM3, RHOt11, RHOt12, RHOt13, ...
    RHOt21, RHOt22, RHOt23, RHOt31, RHOt32, RHOt33];

%% calculate correlations with First Derivative of Spectrum
ll=727; hl=1089;
Y=vis_nir_spectra(1453:1456,1:19)';
X=vis_nir_spectra(ll:hl,1:19)';


for i=1:1:363
        RHOt11(i,1) = corr(X(:,i),Y(:,1));
        RHOt12(i,1) = corr(X(:,i),Y(:,2));
        RHOt13(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,20:38)';

for i=1:1:363
        RHOt21(i,1) = corr(X(:,i),Y(:,1));
        RHOt22(i,1) = corr(X(:,i),Y(:,2));
        RHOt23(i,1) = corr(X(:,i),Y(:,3));
end
  
X=vis_nir_spectra(ll:hl,39:57)';

for i=1:1:363
        RHOt31(i,1) = corr(X(:,i),Y(:,1));
        RHOt32(i,1) = corr(X(:,i),Y(:,2));
        RHOt33(i,1) = corr(X(:,i),Y(:,3));
end

X=avg1(ll:hl,:)';
for i=1:1:363
        RHOM1(i,1) = corr(X(:,i),Y(:,1));
        RHOM2(i,1) = corr(X(:,i),Y(:,2));
        RHOM3(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,:)';
Y=vis_nir_spectra(1453:1456,:)';

for i=1:1:363
        RHOall1(i,1) = corr(X(:,i),Y(:,1));
        RHOall2(i,1) = corr(X(:,i),Y(:,2));
        RHOall3(i,1) = corr(X(:,i),Y(:,3));
end

corr_fd= [RHOall1, RHOall2 RHOall3, RHOM1, RHOM2, RHOM3, RHOt11, RHOt12, RHOt13, ...
    RHOt21, RHOt22, RHOt23, RHOt31, RHOt32, RHOt33];

%% calculate correlations with second derivative Spectrum
ll=1090; hl=1452;
Y=vis_nir_spectra(1453:1456,1:19)';
X=vis_nir_spectra(ll:hl,1:19)';


for i=1:1:363
        RHOt11(i,1) = corr(X(:,i),Y(:,1));
        RHOt12(i,1) = corr(X(:,i),Y(:,2));
        RHOt13(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,20:38)';

for i=1:1:363
        RHOt21(i,1) = corr(X(:,i),Y(:,1));
        RHOt22(i,1) = corr(X(:,i),Y(:,2));
        RHOt23(i,1) = corr(X(:,i),Y(:,3));
end
  
X=vis_nir_spectra(ll:hl,39:57)';

for i=1:1:363
        RHOt31(i,1) = corr(X(:,i),Y(:,1));
        RHOt32(i,1) = corr(X(:,i),Y(:,2));
        RHOt33(i,1) = corr(X(:,i),Y(:,3));
end

X=avg1(ll:hl,:)';
for i=1:1:363
        RHOM1(i,1) = corr(X(:,i),Y(:,1));
        RHOM2(i,1) = corr(X(:,i),Y(:,2));
        RHOM3(i,1) = corr(X(:,i),Y(:,3));
end

X=vis_nir_spectra(ll:hl,:)';
Y=vis_nir_spectra(1453:1456,:)';

for i=1:1:363
        RHOall1(i,1) = corr(X(:,i),Y(:,1));
        RHOall2(i,1) = corr(X(:,i),Y(:,2));
        RHOall3(i,1) = corr(X(:,i),Y(:,3));
end

corr_sd= [RHOall1, RHOall2 RHOall3, RHOM1, RHOM2, RHOM3, RHOt11, RHOt12, RHOt13, ...
    RHOt21, RHOt22, RHOt23, RHOt31, RHOt32, RHOt33];

%% Figures of R2
figure(2)
plot(RSE(1:363,:), corr_raw(:,1).^2, 'or'); hold on
plot(RSE(364:726,:), corr_filt(:,1).^2, 'xb'); hold on
plot(RSE(727:1089,:), corr_fd(:,1).^2, '*k'); hold on
plot(RSE(1090:1452,:), corr_sd(:,1).^2, '+m'); grid off
xlabel('RSE', 'FontSize',16)
ylabel('R ^2, % sand', 'FontSize',16)
set(gca,'XTick',0:1:20); set(gca,'YTick',0.0:0.10:1.0); 
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','NorthEast')
xlim([0 12])
ylim([0.0 0.7])

figure(3)
plot(RSE(1:363,:), corr_raw(:,2).^2, 'or'); hold on
plot(RSE(364:726,:), corr_filt(:,2).^2, 'xb'); hold on
plot(RSE(727:1089,:), corr_fd(:,2).^2, '*k'); hold on
plot(RSE(1090:1452,:), corr_sd(:,2).^2, '+m'); grid off
xlabel('RSE', 'FontSize',16)
ylabel('R ^2, % clay', 'FontSize',16)
set(gca,'XTick',0:1:20); set(gca,'YTick',0.0:0.10:1.0); 
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','NorthEast')
xlim([0 12])
ylim([0.0 0.9])

figure(4)
plot(RSE(1:363,:), corr_raw(:,3).^2, 'or'); hold on
plot(RSE(364:726,:), corr_filt(:,3).^2, 'xb'); hold on
plot(RSE(727:1089,:), corr_fd(:,3).^2, '*k'); hold on
plot(RSE(1090:1452,:), corr_sd(:,3).^2, '+m'); grid off
xlabel('RSE', 'FontSize',16)
ylabel('R ^2, %SOM', 'FontSize',16)
set(gca,'XTick',0:1:20); set(gca,'YTick',0.0:0.10:1.0); 
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','NorthEast')
xlim([0 12])
ylim([0.0 0.7])

figure(5)
plot(wavelengths, avg1(1:363,:), '.-.'); grid off
xlabel('Wavelengths, nm', 'FontSize',16)
ylabel('Reflectance ', 'FontSize',16)
xlim([400 2200])
ylim([0.0 0.6])

%% calculate RMSE predictions with Raw Spectrum

ll=1; hl=363;
Y=vis_nir_spectra(1453:1456,:)';
X=vis_nir_spectra(ll:hl,:)';
c=length(Y);
for i=1:1:363
        

p = polyfit(X(:,i),Y(:,1),1); yhat=p(2)+p(1).*X(:,i); Pall1(i,1)=p(1); RMSEall1(i,1)=sqrt(mean((Y(:,1)-yhat).^2));P2all1(i,1)=p(2);
SDEall1(i,1)=sqrt(sum((Y(:,1)-yhat).^2)/c);% SDE(i,1)=b=yhat-mean(yhat); a= X(:,i)-mean(X(:,i)); sqrt((sum(b.^2) -((sum(a.*b)).^2/(sum(a.^2))))/c);
p = polyfit(X(:,i),Y(:,2),1); yhat=p(2)+p(1).*X(:,i); Pall2(i,1)=p(1); RMSEall2(i,1)=sqrt(mean((Y(:,2)-yhat).^2));P2all2(i,1)=p(2);
SDEall2(i,1)=sqrt(sum((Y(:,2)-yhat).^2)/c);
p = polyfit(X(:,i),Y(:,3),1); yhat=p(2)+p(1).*X(:,i); Pall3(i,1)=p(1); RMSEall3(i,1)=sqrt(mean((Y(:,3)-yhat).^2));P2all3(i,1)=p(2);
SDEall3(i,1)=sqrt(sum((Y(:,3)-yhat).^2)/c);
end

rmse_raw= [RMSEall1, RMSEall2 RMSEall3];
sde_raw=[SDEall1, SDEall2, SDEall3];
slope_raw=[Pall1, Pall2, Pall3];
intcpt_raw=[P2all1, P2all2, P2all3];

%% calculate RMSE predictions with Filtered Spectrum

ll=364; hl=726;
X=vis_nir_spectra(ll:hl,:)';

for i=1:1:363
        
p = polyfit(X(:,i),Y(:,1),1); yhat=p(2)+p(1).*X(:,i); Pall1(i,1)=p(1); RMSEall1(i,1)=sqrt(mean((Y(:,1)-yhat).^2));P2all1(i,1)=p(2);
SDEall1(i,1)=sqrt(sum((Y(:,1)-yhat).^2)/c);
p = polyfit(X(:,i),Y(:,2),1); yhat=p(2)+p(1).*X(:,i); Pall2(i,1)=p(1); RMSEall2(i,1)=sqrt(mean((Y(:,2)-yhat).^2));P2all2(i,1)=p(2);
SDEall2(i,1)=sqrt(sum((Y(:,2)-yhat).^2)/c);
p = polyfit(X(:,i), Y(:,3),1); yhat=p(2)+p(1).*X(:,i); Pall3(i,1)=p(1); RMSEall3(i,1)=sqrt(mean((Y(:,3)-yhat).^2));P2all3(i,1)=p(2);
SDEall3(i,1)=sqrt(sum((Y(:,3)-yhat).^2)/c);
end

rmse_filt= [RMSEall1, RMSEall2 RMSEall3];
sde_filt=[SDEall1, SDEall2, SDEall3];
slope_filt=[Pall1, Pall2, Pall3];
intcpt_filt=[P2all1, P2all2, P2all3];

%% calculate RMSE predictions with First Derivative Spectrum

ll=727; hl=1089;
X=vis_nir_spectra(ll:hl,:)';

for i=1:1:363
        
p = polyfit(X(:,i),Y(:,1),1); yhat=p(2)+p(1).*X(:,i); Pall1(i,1)=p(1); RMSEall1(i,1)=sqrt(mean((Y(:,1)-yhat).^2));P2all1(i,1)=p(2);
SDEall1(i,1)=sqrt(sum((Y(:,1)-yhat).^2)/c);
p = polyfit(X(:,i),Y(:,2),1); yhat=p(2)+p(1).*X(:,i); Pall2(i,1)=p(1); RMSEall2(i,1)=sqrt(mean((Y(:,2)-yhat).^2));P2all2(i,1)=p(2);
SDEall2(i,1)=sqrt(sum((Y(:,2)-yhat).^2)/c);
p = polyfit(X(:,i), Y(:,3),1); yhat=p(2)+p(1).*X(:,i); Pall3(i,1)=p(1); RMSEall3(i,1)=sqrt(mean((Y(:,3)-yhat).^2));P2all3(i,1)=p(2);
SDEall3(i,1)=sqrt(sum((Y(:,3)-yhat).^2)/c);
end


rmse_fd= [RMSEall1, RMSEall2 RMSEall3];
sde_fd=[SDEall1, SDEall2, SDEall3];
slope_fd=[Pall1, Pall2, Pall3];
intcpt_fd=[P2all1, P2all2, P2all3];

%% calculate RMSE predictions with Filtered Spectrum

ll=1090; hl=1452;
X=vis_nir_spectra(ll:hl,:)';

for i=1:1:363
        
p = polyfit(X(:,i),Y(:,1),1); yhat=p(2)+p(1).*X(:,i); Pall1(i,1)=p(1); RMSEall1(i,1)=sqrt(mean((Y(:,1)-yhat).^2));P2all1(i,1)=p(2);
SDEall1(i,1)=sqrt(sum((Y(:,1)-yhat).^2)/c);
p = polyfit(X(:,i),Y(:,2),1); yhat=p(2)+p(1).*X(:,i); Pall2(i,1)=p(1); RMSEall2(i,1)=sqrt(mean((Y(:,2)-yhat).^2));P2all2(i,1)=p(2);
SDEall2(i,1)=sqrt(sum((Y(:,2)-yhat).^2)/c);
p = polyfit(X(:,i), Y(:,3),1); yhat=p(2)+p(1).*X(:,i); Pall3(i,1)=p(1); RMSEall3(i,1)=sqrt(mean((Y(:,3)-yhat).^2));P2all3(i,1)=p(2);
SDEall3(i,1)=sqrt(sum((Y(:,3)-yhat).^2)/c);
end

rmse_sd= [RMSEall1, RMSEall2 RMSEall3];
sde_sd=[SDEall1, SDEall2, SDEall3];
slope_sd=[Pall1, Pall2, Pall3];
intcpt_sd=[P2all1, P2all2, P2all3];


figure(6)
plot(abs(slope_raw(:,1)).*RMSE(1:363,:),rmse_raw(:,1), 'or'); hold on
plot(abs(slope_filt(:,1)).*RMSE(364:726,:),rmse_filt(:,1), 'xb'); hold on
plot(abs(slope_fd(:,1)).*RMSE(727:1089,:),rmse_fd(:,1), '*k'); hold on
plot(abs(slope_sd(:,1)).*RMSE(1090:1452,:),rmse_sd(:,1), '+m'); hold on
zz=std(Y(:,1)); plot((0:2:10), [zz,zz,zz,zz,zz,zz], 'linewidth',3); grid off
xlabel('Measurement Precision, % sand', 'FontSize',16)
ylabel('Standard Error of Prediction, % sand', 'FontSize',16)
set(gca,'XTick',0:1:20); set(gca,'YTick',0.0:2:20); 
plot([15,0], [0,15], '.-.'); plot([16,0], [0,16], '.-.')
plot([17,0], [0,17], '.-.'); plot([18,0], [0,18], '.-.')
plot([19,0], [0,19], '.-.'); plot([20,0], [0,20], '.-.')
plot([21,0], [0,21], '.-.');plot([22,0], [0,22], '.-.')
plot([23,0], [0,23], '.-.')
plot(3.12, 4.82, 'dg', 'markersize',12)
plot(2.64, 4.09, 'or', 'markersize',12)
plot(3.43, 4.96, 'xb', 'markersize',12)
plot(3.29, 4.46, '*k', 'markersize',12)
plot(3.09, 5.85, '+m', 'markersize',12)
xlim([0 10]); ylim([0 20])
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','SouthEast')

figure(7)
plot(abs(slope_raw(:,2)).*RMSE(1:363,:),rmse_raw(:,2), 'or'); hold on
plot(abs(slope_filt(:,2)).*RMSE(364:726,:),rmse_filt(:,2), 'xb'); hold on
plot(abs(slope_fd(:,2)).*RMSE(727:1089,:),rmse_fd(:,2), '*k'); hold on
plot(abs(slope_sd(:,2)).*RMSE(1090:1452,:),rmse_sd(:,2), '+m'); hold on
zz=std(Y(:,2)); plot((0:5), [zz,zz,zz,zz,zz zz], 'linewidth',3); grid off
xlabel('Measurement Precision, % clay', 'FontSize',16)
ylabel('Standard Error of Prediction, % clay', 'FontSize',16)
set(gca,'XTick',0:1:10); set(gca,'YTick',0:3:20); 
plot([6,0], [0,6], '.-.'); plot([7,0], [0,7], '.-.')
plot([8,0], [0,8], '.-.'); plot([9,0], [0,9], '.-.')
plot([10,0], [0,10], '.-.'); plot([11,0], [0,11], '.-.')
plot([12,0], [0,12], '.-.'); plot([5,0], [0,5], '.-.');
plot(1.81, 3.48, 'dg', 'markersize',12)
plot(1.91, 4.59, 'or', 'markersize',12)
plot(1.66, 2.89, 'xb', 'markersize',12)
plot(1.28, 1.73, '*k', 'markersize',12)
plot(1.23, 2.04, '+m', 'markersize',12)
xlim([0 5]); ylim([0 12])
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','SouthEast')

figure(8)
plot(abs(slope_raw(:,3)).*RMSE(1:363,:),rmse_raw(:,3), 'or'); hold on
plot(abs(slope_filt(:,3)).*RMSE(364:726,:),rmse_filt(:,3), 'xb'); hold on
plot(abs(slope_fd(:,3)).*RMSE(727:1089,:),rmse_fd(:,3), '*k'); hold on
plot(abs(slope_sd(:,3)).*RMSE(1090:1452,:),rmse_sd(:,3), '+m'); hold on
zz=std(Y(:,3)); plot((0:4), [zz,zz,zz,zz,zz], 'linewidth',3); grid off
xlabel('Measurement Precision, % SOM', 'FontSize',16)
ylabel('Standard Error of Prediction, % SOM', 'FontSize',16)
set(gca,'XTick',0:1:10); set(gca,'YTick',0:3:20); 
plot([4,0], [0,4], '.-.');plot([4.5,0], [0,4.5], '.-.')
plot([5,0], [0,5], '.-.') ;plot([5.5,0], [0,5.5], '.-.'); 
plot([6,0], [0,6], '.-.')
plot(1.38, 2.14, 'dg', 'markersize',12)
plot(1.62, 4.18, 'or', 'markersize',12)
plot(1.67, 2.44, 'xb', 'markersize',12)
plot(1.31, 2.05, '*k', 'markersize',12)
plot(1.39, 2.75, '+m', 'markersize',12)
xlim([0 3]); ylim([0 7])
legend('Raw',  'Smooth', '1^s^t SGD', '2^n^d SGD', ...
    'Location','SouthEast')

%[minima1.DataIndex, wavelengths(minima1.DataIndex), minima1.Position(2)+ minima1.Position(1), minima1.Position(2), minima1.Position(2)- minima1.Position(1)]

figure(9)
tt= mean(avg1(1:363,:),2);
tt1=std(avg1(1:363,:)')';
plot(wavelengths, tt+tt1, '--g', 'markersize',5); hold on
plot(wavelengths, tt, '-k', 'markersize',5); hold on
plot(wavelengths, tt-tt1, '--g', 'markersize',5); hold on

text(1462, tt(190),'\leftarrow sand (best)','HorizontalAlignment','left')

text(1940,tt(296),' sand, clay (best) \rightarrow', 'HorizontalAlignment','right')

text(421,tt(2),'\leftarrow SOM ', 'HorizontalAlignment','left')
text(1747,tt(251),'\leftarrow SOM (best) ', 'HorizontalAlignment','left')
text(1605,tt(220),'SOM \rightarrow  ', 'HorizontalAlignment','right')
text(590,tt(31),'\leftarrow SOM ', 'HorizontalAlignment','left')
text(722,tt(55),'\leftarrow SOM ', 'HorizontalAlignment','left')

set(gca,'XTick',400:300:2200); set(gca,'YTick',0:0.1:0.6); 
xlabel('Wavelengths, nm', 'FontSize',16)
ylabel('Reflectance ', 'FontSize',16)
xlim([400 2200])
ylim([0.0 0.52])
grid off


%% figure (6);
% ll=1;hl=363;
% 
% subplot(2,2,1);plot(vis_nir_factors(ll:hl,:), RSE(ll:hl,:), '.m');
% xlim([min(wavelengths) max(wavelengths)])
% ylim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('RSE')
% grid off
% 
% subplot(2,2,2); plot(wavelengths,corr_raw(:,1),  '.r');hold on
% subplot(2,2,2);plot(wavelengths,corr_raw(:,4),  '+g');hold on
% subplot(2,2,2); plot(wavelengths,corr_raw(:,7),  'xb');hold on
% subplot(2,2,2);plot(wavelengths,corr_raw(:,10),  '*m');hold on
% subplot(2,2,2);plot(wavelengths,corr_raw(:,13),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-SAND');%zlabel('RSE')
% title('Raw VIS-NIR Spectrum, correlation coefficients (R) and RSE')
% % legend ALL MEAN T1 T2 T3
% 
% subplot(2,2,3); plot(wavelengths,corr_raw(:,2),  '.r');hold on
% subplot(2,2,3);plot(wavelengths,corr_raw(:,5),  '+g');hold on
% subplot(2,2,3); plot(wavelengths,corr_raw(:,8),  'xb');hold on
% subplot(2,2,3);plot(wavelengths,corr_raw(:,11),  '*m');hold on
% subplot(2,2,3);plot(wavelengths,corr_raw(:,14),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CLAY');%zlabel('RSE')
% 
% 
% subplot(2,2,4); plot(wavelengths,corr_raw(:,3),  '.r');hold on
% subplot(2,2,4);plot(wavelengths,corr_raw(:,6),  '+g');hold on
% subplot(2,2,4); plot(wavelengths,corr_raw(:,9),  'xb');hold on
% subplot(2,2,4);plot(wavelengths,corr_raw(:,12),  '*m');hold on
% subplot(2,2,4);plot(wavelengths,corr_raw(:,15),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CARBON');%zlabel('RSE')
% 
% 
% figure (7);
% ll=364;hl=726;
% 
% subplot(2,2,1);plot(vis_nir_factors(ll:hl,:), RSE(ll:hl,:), '.b');
% xlim([min(wavelengths) max(wavelengths)])
% ylim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('RSE')
% grid off
% 
% subplot(2,2,2); plot(wavelengths,corr_filt(:,1),  '.r');hold on
% subplot(2,2,2);plot(wavelengths,corr_filt(:,4),  '+g');hold on
% subplot(2,2,2); plot(wavelengths,corr_filt(:,7),  'xb');hold on
% subplot(2,2,2);plot(wavelengths,corr_filt(:,10), '*m');hold on
% subplot(2,2,2);plot(wavelengths,corr_filt(:,13),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-SAND');%zlabel('RSE')
% title('Denoised VIS-NIR Spectrum, correlation coefficients (R) and RSE')
% % legend ALL MEAN T1 T2 T3
% 
% subplot(2,2,3); plot(wavelengths,corr_filt(:,2),  '.r');hold on
% subplot(2,2,3);plot(wavelengths,corr_filt(:,5),  '+g');hold on
% subplot(2,2,3); plot(wavelengths,corr_filt(:,8),  'xb');hold on
% subplot(2,2,3);plot(wavelengths,corr_filt(:,11),  '*m');hold on
% subplot(2,2,3);plot(wavelengths,corr_filt(:,14),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CLAY');%zlabel('RSE')
% 
% 
% subplot(2,2,4); plot(wavelengths,corr_filt(:,3),  '.r');hold on
% subplot(2,2,4);plot(wavelengths,corr_filt(:,6), '+g');hold on
% subplot(2,2,4); plot(wavelengths,corr_filt(:,9),  'xb');hold on
% subplot(2,2,4);plot(wavelengths,corr_filt(:,12),  '*m');hold on
% subplot(2,2,4);plot(wavelengths,corr_filt(:,15),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CARBON');%zlabel('RSE')
% 
% figure (8);
% ll=727;hl=1089;
% 
% subplot(2,2,1);plot(vis_nir_factors(ll:hl,:), RSE(ll:hl,:), '.g');
% xlim([min(wavelengths) max(wavelengths)])
% ylim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('RSE')
% grid off
% 
% subplot(2,2,2); plot(wavelengths,corr_fd(:,1),  '.r');hold on
% subplot(2,2,2);plot(wavelengths,corr_fd(:,4),  '+g');hold on
% subplot(2,2,2); plot(wavelengths,corr_fd(:,7),  'xb');hold on
% subplot(2,2,2);plot(wavelengths,corr_fd(:,10),  '*m');hold on
% subplot(2,2,2);plot(wavelengths,corr_fd(:,13),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-SAND');%zlabel('RSE')
% title('Ist Derivative of VIS-NIR Spectrum, correlation coefficients (R) and RSE')
% % legend ALL MEAN T1 T2 T3
% 
% subplot(2,2,3); plot(wavelengths,corr_fd(:,2), '.r');hold on
% subplot(2,2,3);plot(wavelengths,corr_fd(:,5),  '+g');hold on
% subplot(2,2,3); plot(wavelengths,corr_fd(:,8),  'xb');hold on
% subplot(2,2,3);plot(wavelengths,corr_fd(:,11),  '*m');hold on
% subplot(2,2,3);plot(wavelengths,corr_fd(:,14),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CLAY');%zlabel('RSE')
% 
% 
% subplot(2,2,4); plot(wavelengths,corr_fd(:,3),  '.r');hold on
% subplot(2,2,4);plot(wavelengths,corr_fd(:,6),  '+g');hold on
% subplot(2,2,4); plot(wavelengths,corr_fd(:,9),  'xb');hold on
% subplot(2,2,4);plot(wavelengths,corr_fd(:,12),  '*m');hold on
% subplot(2,2,4);plot(wavelengths,corr_fd(:,15),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CARBON');%zlabel('RSE')
% 
% 
% figure (9);
% ll=1090;hl=1452;
% 
% subplot(2,2,1);plot(vis_nir_factors(ll:hl,:), RSE(ll:hl,:), '.r');
% xlim([min(wavelengths) max(wavelengths)])
% ylim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('RSE')
% grid off
% 
% subplot(2,2,2); plot(wavelengths,corr_sd(:,1),  '.r');hold on
% subplot(2,2,2);plot(wavelengths,corr_sd(:,4), '+g');hold on
% subplot(2,2,2); plot(wavelengths,corr_sd(:,7), 'xb');hold on
% subplot(2,2,2);plot(wavelengths,corr_sd(:,10), '*m');hold on
% subplot(2,2,2);plot(wavelengths,corr_sd(:,13), 'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-SAND');%zlabel('RSE')
% title('IInd Derivative VIS-NIR Spectrum, correlation coefficients (R) and RSE')
% % legend ALL MEAN T1 T2 T3
% 
% subplot(2,2,3); plot(wavelengths,corr_sd(:,2), '.r');hold on
% subplot(2,2,3);plot(wavelengths,corr_sd(:,5),  '+g');hold on
% subplot(2,2,3); plot(wavelengths,corr_sd(:,8),  'xb');hold on
% subplot(2,2,3);plot(wavelengths,corr_sd(:,11),  '*m');hold on
% subplot(2,2,3);plot(wavelengths,corr_sd(:,14),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CLAY');%zlabel('RSE')
% 
% 
% subplot(2,2,4); plot(wavelengths,corr_sd(:,3), '.r');hold on
% subplot(2,2,4);plot(wavelengths,corr_sd(:,6),  '+g');hold on
% subplot(2,2,4); plot(wavelengths,corr_sd(:,9),  'xb');hold on
% subplot(2,2,4);plot(wavelengths,corr_sd(:,12),  '*m');hold on
% subplot(2,2,4);plot(wavelengths,corr_sd(:,15),  'ok');hold on
% grid off
% xlim([min(wavelengths) max(wavelengths)])
% ylim([-1.0 1.0])
% %zlim([min(RSE(ll:hl,:)) max(RSE(ll:hl,:))])
% xlabel('Wavelengths, nm'); ylabel('correlations(R)-CARBON');%zlabel('RSE')