% MATLAB code for generating Figure 4 of the manuscript
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
        [L0a,L0b,L1a,L1b,L0L0a,L0L0b,L0L1a,L1L0a,L1L1a,L1L0b,L0L1b,L1L1b] = Kolmoperator(aT,bT,af,bf,zvar);
        L0Af(i,j) = L0a;
        L0Bf(i,j) = L0b;
        L1Af(i,j) = L1a;
        L1Bf(i,j) = L1b;
        L0L0Af(i,j) = L0L0a;
        L0L0Bf(i,j) = L0L0b;
        L0L1Af(i,j) = L0L1a;
        L1L0Af(i,j) = L1L0a;
        L1L1Af(i,j) = L1L1a;
        L1L0Bf(i,j) = L1L0b;
        L0L1Bf(i,j) = L0L1b;
        L1L1Bf(i,j) = L1L1b;
        clear aT bT
    end
end
%
xi = [z1; z2; z3];
av = matlabFunction(af,'vars',{([xi])});
bv = matlabFunction(bf,'vars',{([xi])});
L0Av = matlabFunction(L0Af,'vars',{([xi])});
L1Av = matlabFunction(L1Af,'vars',{([xi])});
L0Bv = matlabFunction(L0Bf,'vars',{([xi])});
L1Bv = matlabFunction(L1Bf,'vars',{([xi])});
L1L1Bv = matlabFunction(L1L1Bf,'vars',{([xi])});
L1L0Bv = matlabFunction(L1L0Bf,'vars',{([xi])});
L0L0Bv = matlabFunction(L0L0Bf,'vars',{([xi])});
L0L1Bv = matlabFunction(L0L1Bf,'vars',{([xi])});
L1L1Av = matlabFunction(L1L1Af,'vars',{([xi])});
L1L0Av = matlabFunction(L1L0Af,'vars',{([xi])});
L0L0Av = matlabFunction(L0L0Af,'vars',{([xi])});
L0L1Av = matlabFunction(L0L1Af,'vars',{([xi])});
%
% Formulation considering Euclidean setting
aEf = az*zvar;
bEf = bf*zvar;
%
[L0AEf,L0BEf,L1AEf,L1BEf,L0L0AEf,L0L0BEf,L0L1AEf,L1L0AEf,L1L1AEf,L1L0BEf,L0L1BEf,L1L1BEf] = tayL0L1_w3(aEf,bEf,zvar);
%
aEv = matlabFunction(aEf,'vars',{([xi])});
bEv = matlabFunction(bEf,'vars',{([xi])});
L0AEv = matlabFunction(L0AEf,'vars',{([xi])});
L1AEv = matlabFunction(L1AEf,'vars',{([xi])});
L0BEv = matlabFunction(L0BEf,'vars',{([xi])});
L1BEv = matlabFunction(L1BEf,'vars',{([xi])});
L1L1BEv = matlabFunction(L1L1BEf,'vars',{([xi])});
L1L0BEv = matlabFunction(L1L0BEf,'vars',{([xi])});
L0L0BEv = matlabFunction(L0L0BEf,'vars',{([xi])});
L0L1BEv = matlabFunction(L0L1BEf,'vars',{([xi])});
L1L1AEv = matlabFunction(L1L1AEf,'vars',{([xi])});
L1L0AEv = matlabFunction(L1L0AEf,'vars',{([xi])});
L0L0AEv = matlabFunction(L0L0AEf,'vars',{([xi])});
L0L1AEv = matlabFunction(L0L1AEf,'vars',{([xi])});
%
dt = 0.01; % Time discretization
t = 0:dt:T;
for scheme = 1:6

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
for l=2:numel(t)    
  dB(l) = sqrt(dt)*randn;
  %
    dB1(l)=sqrt(dt)*randn;
    dB2(l)=sqrt(dt)*randn;
% 
% MSIs
    Iw=dB(l);
    Iws=0.5*dt*(dB1(l)+dB2(l)/sqrt(3));
    Isw=0.5*dt*(dB1(l)-dB2(l)/sqrt(3));
    Iww=0.5*(Iw*Iw-dt);
    Iss = 0.5*(dt)^2;
    Iwww=(1/6)*(Iw*Iw-3*dt)*Iw;
    Iwss=1/6*Iw*dt^2;
    Isww = 1/6*(Iw*Iw-dt)*dt;
    Isss=1/6*dt^3;
%   
    a = av([z(:,l-1)]);
    b = bv([z(:,l-1)]);
    L0A = L0Av([z(:,l-1)]);
    L1A = L1Av([z(:,l-1)]);
    L0B = L0Bv([z(:,l-1)]);
    L1B = L1Bv([z(:,l-1)]);
    L1L1B = L1L1Bv([z(:,l-1)]);
    L1L0B = L1L0Bv([z(:,l-1)]);
    L0L1B = L0L1Bv([z(:,l-1)]);
    L0L0B = L0L0Bv([z(:,l-1)]);
    L1L1A = L1L1Av([z(:,l-1)]);
    L1L0A = L1L0Av([z(:,l-1)]);
    L0L1A = L0L1Av([z(:,l-1)]);
    L0L0A = L0L0Av([z(:,l-1)]);
%
    aE = aEv([z(:,l-1)]);
    bE = bEv([z(:,l-1)]);
    L0EA = L0AEv([z(:,l-1)]);
    L1EA = L1AEv([z(:,l-1)]);
    L0EB = L0BEv([z(:,l-1)]);
    L1EB = L1BEv([z(:,l-1)]);
    L1L1EB = L1L1BEv([z(:,l-1)]);
    L1L0EB = L1L0BEv([z(:,l-1)]);
    L0L1EB = L0L1BEv([z(:,l-1)]);
    L0L0EB = L0L0BEv([z(:,l-1)]);
    L1L1EA = L1L1AEv([z(:,l-1)]);
    L1L0EA = L1L0AEv([z(:,l-1)]);
    L0L1EA = L0L1AEv([z(:,l-1)]);
    L0L0EA = L0L0AEv([z(:,l-1)]);
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
        for ii = 1:length(a)
            for jj = 1:length(a)
             Omega(ii,jj) = Omegzero(ii,jj) + a(ii,jj)*dt + b(ii,jj)*Iw + L0A(ii,jj)*Iss + L1B(ii,jj)*Iww...
                             + L1A(ii,jj)*Iws + L0B(ii,jj)*Isw; % g-Weak 2.0
             end
        end
         z(:,l) = so3_exp_new(Omega)*z(:,l-1);
%
elseif scheme == 3    
        for ii = 1:length(a)
            for jj = 1:length(a)
             Omega(ii,jj) = Omegzero(ii,jj) + a(ii,jj)*dt + b(ii,jj)*Iw + L0A(ii,jj)*Iss + L1B(ii,jj)*Iww...
                             + L1A(ii,jj)*Iws + L0B(ii,jj)*Isw + L0L0A(ii,jj)*Isss + (L0L0B(ii,jj) + L0L1A(ii,jj) + L1L0A(ii,jj))*Iwss...
                             + (L1L1A(ii,jj) + L1L0B(ii,jj) + L0L1B(ii,jj))*Isww + L1L1B(ii,jj)*Iwww; % g-Weak 3.0
             end
         end
        z(:,l) = so3_exp_new(Omega)*z(:,l-1);
%
elseif scheme == 4
        z(:,l) = z(:,l-1) + aE*dt + bE*Iw; % Weak 1.0
%
elseif scheme == 5
        z(:,l) = z(:,l-1) + aE*dt + bE*Iw + L0EA*Iss + L1EB*Iww...
                             + L1EA*Iws + L0EB*Isw; % Weak 2.0
%
elseif scheme == 6
        z(:,l) = z(:,l-1) + aE*dt + bE*Iw + L0EA*Iss + L1EB*Iww...
                             + L1EA*Iws + L0EB*Isw + L0L0EA*Isss + (L0L0EB + L0L1EA + L1L0EA)*Iwss...
                             + (L1L1EA + L1L0EB + L0L1EB)*Isww + L1L1EB*Iwww; % Weak 3.0
end
%
end
%
Q(:,k,:) = z;
clear z
end
    if scheme == 1
        for ii = 1:length(Q)
            Y_GW1C(:,ii) = karcher_mean_sphere(Q(:,:,ii));
        end
    elseif scheme == 2
        for ii = 1:length(Q)
            Y_GW2C(:,ii) = karcher_mean_sphere(Q(:,:,ii));
        end
    elseif scheme == 3
        for ii = 1:length(Q)
            Y_GW3C(:,ii) = karcher_mean_sphere(Q(:,:,ii));
        end
    elseif scheme == 4
        for ii = 1:length(Q)
            Y_GW1IC(:,ii) = mean(Q(:,:,ii),2);
        end
    elseif scheme == 5
        for ii = 1:length(Q)
            Y_GW2IC(:,ii) = mean(Q(:,:,ii),2);
        end
    elseif scheme == 6
        for ii = 1:length(Q)
            Y_GW3IC(:,ii) = mean(Q(:,:,ii),2);
        end
%
    end
clear Q
end
%
%% Plotting
Q_dot1 = abs(dot(Y_GW1C,Y_GW1C)-1);
Q_dot2 = abs(dot(Y_GW2C,Y_GW2C)-1);
Q_dot3 = abs(dot(Y_GW3C,Y_GW3C)-1);
Q_dotI1 = abs(dot(Y_GW1IC,Y_GW1IC)-1);
Q_dotI2 = abs(dot(Y_GW2IC,Y_GW2IC)-1);
Q_dotI3 = abs(dot(Y_GW3IC,Y_GW3IC)-1);
for i = 5:length(Y_GW1C)
if Q_dot1(i) == 0
    Q_dot1(i) = Q_dot1(i-1);
end
if Q_dot2(i) == 0
    Q_dot2(i) = Q_dot2(i-1);
end
if Q_dot3(i) == 0
    Q_dot3(i) = Q_dot3(i-1);
end
end
%
% Create figure
figure1 = figure;
%
% Create axes
axes1 = axes('Position',...
    [0.106382978723404 0.120325811399254 0.842062193126023 0.804674188600746]);
hold(axes1,'on');
%
figure(1);
semilogy((Q_dotI1(2:end)),'linewidth',2,'LineStyle','--','DisplayName','Weak 1.0'); hold on
semilogy((Q_dotI2(2:end)),'linewidth',2,'LineStyle','--','DisplayName','Weak 2.0'); hold on
semilogy((Q_dotI3(2:end)),'linewidth',2,'LineStyle','--','DisplayName','Weak 3.0'); hold on
semilogy((Q_dot1(2:end)),'linewidth',2,'DisplayName','g-Weak 1.0'); hold on
semilogy((Q_dot2(2:end)),'linewidth',2,'DisplayName','g-Weak 2.0'); hold on
semilogy((Q_dot3(2:end)),'linewidth',2,'DisplayName','g-Weak 3.0'); hold on
hold off
ylabel('$|1-\|z\||$','FontWeight','bold','FontSize',30,...
    'FontName','Times New Roman',...
    'Interpreter','latex')
xlabel('Samples','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
ylim([10^(-18) 10])
xlim([0 200])
%
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',22,'FontWeight','bold',...
    'XTick',[0 40 80 120 160 200],'XTickLabel',...
    {'0','40','80','120','160','200'},'YScale','log');
%
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.826067440426325 0.41291202719847 0.098405715330018 0.228622321207563]);
%
% Create textbox
annotation(figure1,'textbox',...
    [0.108019639934533 0.173930907110063 0.839607201309328 0.123068320158762],...
    'FitBoxToText','off',...
    'FaceAlpha',0.2,...
    'EdgeColor','none',...
    'BackgroundColor',[0.650980392156863 0.650980392156863 0.650980392156863]);
%
% Create textarrow
annotation(figure1,'textarrow',[0.400109769484083 0.445584744928622],...
    [0.432703003337041 0.30366179781583],'String',{'Geometric schemes'},...
    'LineWidth',2,...
    'HeadStyle','cback2',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Times New Roman');
%
% Create textarrow
annotation(figure1,'textarrow',[0.680851063829787 0.746934940551963],...
    [0.596482233315753 0.701369451360867],'String',{'Non-geometric schemes'},...
    'LineWidth',2,...
    'HeadStyle','cback2',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Times New Roman');
%
% Create textbox
annotation(figure1,'textbox',...
    [0.107905039096738 0.701278594609927 0.840425531914893 0.22],...
    'FitBoxToText','off',...
    'FaceAlpha',0.15,...
    'EdgeColor','none',...
    'BackgroundColor',[0 1 1]);
%
% End