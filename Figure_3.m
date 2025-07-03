% MATLAB code for generating Figure 3 of the manuscript
clc
clear
close all
warning off
%
T = 2;    % total time of integration
%
N=10000;   % samples in each MC run
%
% Parameters of Stochastic Rigid Body
I1 = 2; I2 = 0.4; I3 = 0.6; %
zzero = [sin(1.1); 0; cos(1.1)];   % initial condition
syms z1 z2 z3
%
zvar = [z1; z2; z3];
sigmaz = [0          -z3/I3      z2/I2;
          z3/I3      0          -z1/I1;
         -z2/I2   z1/I1          0];
%
az = [-1/2*z2^2/I2^2-1/2*z3^2/I3^2                              0                                                      0;
                    (z1*z2)/(I1*I2)                        -1/2*z1^2/I1^2-1/2*z3^2/I3^2                                0;
                    (z1*z3)/(I1*I3)                                        (z2*z3)/(I2*I3)                 -1/2*z1^2/I1^2-1/2*z2^2/I2^2];
af = az-1/2*sigmaz^2;
bf = sigmaz;
%
% Formulation for geometric SDE
for i = 1:length(af)
    for j = 1:length(af)
        aT = af(i,j);
        bT = bf(i,j);
        clear aT bT
    end
end
%
xi = [z1; z2; z3];
av = matlabFunction(af,'vars',{([xi])});
bv = matlabFunction(bf,'vars',{([xi])});
%
% Formulation considering Euclidean setting
aEf = az*zvar;
bEf = bf*zvar;
%
aEv = matlabFunction(aEf,'vars',{([xi])});
bEv = matlabFunction(bEf,'vars',{([xi])});
%
dt = 0.01; % Time discretization
t = 0:dt:T;
%
for scheme = 1:2
%
for k = 1:N
    k
%
% On the tangent space
Omegzero = zeros(length(af),length(af));
Omega = zeros(length(af),length(af));
%
z = zeros(3,numel(t));
z(:,1) = zzero;
%
rng(0);
for l=2:numel(t)    
  dB(l) = sqrt(dt)*randn;
% Weak 3.0 solution
 %
    dB1(l)=sqrt(dt)*randn;
    dB2(l)=sqrt(dt)*randn;
% 
% MSIs
    Iw=dB(l); 
%   
    a = av([z(:,l-1)]);
    b = bv([z(:,l-1)]);
%
    aE = aEv([z(:,l-1)]);
    bE = bEv([z(:,l-1)]);
 %
    if scheme == 1
        for ii = 1:length(a)
            for jj = 1:length(a)
             Omega(ii,jj) = Omegzero(ii,jj) + a(ii,jj)*dt + b(ii,jj)*Iw; % g-Weak 1.0
            end
         end
        z(:,l) = so3_exp_new(Omega)*z(:,l-1);
 %
    elseif scheme == 2
        z(:,l) = z(:,l-1) + aE*dt + bE*Iw; % Weak 1.0
    end
end
%
Q(:,k,:) = z;
clear z
end
%
    if scheme == 1
        for ii = 1:length(Q)
        Y_GW1C(:,ii) = karcher_mean_sphere(Q(:,:,ii));
        end
    elseif scheme == 2
        for ii = 1:length(Q)
        Y_GW1IC(:,ii) = karcher_mean_sphere(Q(:,:,ii));
        end
    end
%
end
%
%% Plotting
[X1,Y1,Z1] = sphere(25);
X2 = X1;
Y2 = Y1;
Z2 = Z1;
figure(1),
%
% Create axes
axes1 = axes('Position',...
    [0.112252007853986 0.101503927016137 0.771119126263787 0.914295226775633]);
hold(axes1,'on');
%
FIG = surf(X2,Y2,Z2,'DisplayName','S^2','FaceAlpha',0.5,...
    'EdgeAlpha',0.2);
set(FIG, 'FaceAlpha', 0.5,'EdgeAlpha', 0.2)
% shading interp
axis equal
hold on
plot3(Y_GW1C(1,2:end), Y_GW1C(2,2:end), Y_GW1C(3,2:end),'linewidth',3); hold on
plot3(Y_GW1IC(1,2:end), Y_GW1IC(2,2:end), Y_GW1IC(3,2:end),'LineStyle','--','linewidth',3)
hold on
view(1.063829218171949e+02,10.953605815995953)
xlim([-0.2 1])
ylim([-1 0.4])
zlim([0 1])
% Create zlabel
zlabel({'Z-position'});
% Create ylabel
ylabel('Y-position');
% Create xlabel
xlabel('X-position','Rotation',45,'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
%
view(axes1,[106.382921817195 10.953605815996]);
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'DataAspectRatio',[1 1 1],'FontName','Times New Roman','FontSize',...
    30,'FontWeight','bold','LineWidth',1,'XLimitMethod','tight','XTick',...
    [-0.2 0.2 0.6 1],'YLimitMethod','tight','YTick',[-1 -0.6 -0.2 0.2 0.6],...
    'ZLimitMethod','tight');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.111929435314602 0.668134686886139 0.211206891567543 0.194904452979944],...
    'FontSize',30);
legend("S^2","g-Weak 1.0","Weak 1.0")
%
% End