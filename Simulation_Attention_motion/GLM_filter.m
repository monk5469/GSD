function y_filter = GLM_filter (y_noise,DCM)

%%%%%%%%%% GLM model %%%%%%%%%%%%
nu      = size (DCM.U.u,2);
H       = repmat (fft(DCM.h),1,nu);
r_dt    = DCM.Y.dt/DCM.U.dt;
utmp    = full(DCM.U.u)/r_dt;
UH      = fft(utmp).*H;
uh      = ifft(UH);

uh_down = uh(1:r_dt:end,:);

Ny     = size(y_noise,1);
Inter  = 0.5;
NCos   = 180;
SCos = cos(2*pi*(0:1:Ny-1)'*(0:Inter:NCos*Inter)/Ny);
Beta   = pinv([uh_down SCos])*y_noise;


y_filter = [uh_down SCos]*Beta;
% y_test   = [uh_down ]*Beta(1:50,:);

% NoR = 3;
% figure(1);
% plot(y_noise(:,NoR),'g');
% hold on;
% plot(y_filter(:,NoR),'r');
% hold on;
% plot(DCM.y(:,NoR),'b');

% X_design = [y_filter(1:end-1,:) uh_down(1:end-1,:)];
% est       = pinv(X_design)*(diff(y_filter)/DCM.Y.dt);
% 
% figure(2);
% % plot(abs(fft(y_noise(:,NoR))),'g');
% % hold on;
% plot(abs(fft(y_filter(:,NoR))),'r');
% hold on;
% plot(abs(fft(DCM.y(:,NoR))),'b');



