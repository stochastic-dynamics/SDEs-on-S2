% MATLAB code for generating Figure 5 of the manuscript
clc
clear
close all
warning off
%
%% 
syms z [3,1]; 
T = 1;    % total time of integration
%
rng(25); 
A = -eye(3) + ones(3);%
bf =  so3_wedge([0.002; 0.002; 0.002]);
af = (so3_wedge(ones(3,1)) + (1-z.'*z)*A - so3_wedge(cross(A*z,z)));
zzero = [0; 1; 0];   % initial condition
%
% Formulation for geometric SDE
for i = 1:length(af)
%
    for j = 1:length(af)
        aT = af(i,j);
        bT = bf(i,j);
        [L0a,L0b,L1a,L1b,L0L0a,L0L0b,L0L1a,L1L0a,L1L1a,L1L0b,L0L1b,L1L1b] = Kolmoperator(aT,bT,af,bf,z);
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
bv = bf;
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
for int=1:11  % integration with different time steps
%
int
  if(int==1); dt=2^-16;
    elseif(int==2); dt=2^-12;
    elseif(int==3); dt=2^-11;
    elseif(int==4); dt=2^-10;
    elseif(int==5); dt=2^-9;
    elseif(int==6); dt=2^-8;
    elseif(int==7); dt=2^-7;
    elseif(int==8); dt=2^-6;
    elseif(int==9); dt=2^-5;
    elseif(int==10); dt=2^-4;
    elseif(int==11); dt=2^-3;
  end
%
    t = 0:dt:T;
    if int == 1
        SCHEME = 3;
    else
        SCHEME = 1:3;
    end
for scheme = SCHEME
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
    b = bf;
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
for ii = 1:length(a)
    for jj = 1:length(a)
%
    if scheme == 1
        Omega(ii,jj) = Omegzero(ii,jj) + a(ii,jj)*dt + b(ii,jj)*Iw + L0A(ii,jj)*Iss + L1B(ii,jj)*Iww...
                                     + L1A(ii,jj)*Iws + L0B(ii,jj)*Isw + L0L0A(ii,jj)*Isss + (L0L0B(ii,jj) + L0L1A(ii,jj) + L1L0A(ii,jj))*Iwss...
                                     + (L1L1A(ii,jj) + L1L0B(ii,jj) + L0L1B(ii,jj))*Isww + L1L1B(ii,jj)*Iwww; % g-Weak 3.0
    elseif scheme == 2
        Omega(ii,jj) = Omegzero(ii,jj) + a(ii,jj)*dt + b(ii,jj)*Iw + L0A(ii,jj)*Iss + L1B(ii,jj)*Iww...
                                     + L1A(ii,jj)*Iws + L0B(ii,jj)*Isw; % g-Weak 2.0
    elseif scheme == 3
        Omega(ii,jj) = Omegzero(ii,jj) + a(ii,jj)*dt + b(ii,jj)*Iw; % g-Weak 1.0
    end
    end
end
%
z(:,l) = so3_exp_new(Omega)*z(:,l-1);
%
end
%
Q = z;
clear z
%
    if scheme == 1
        Y_GW3 = Q(:,end);
    elseif scheme == 2
        Y_GW2 = Q(:,end);
    elseif scheme == 3
        Y_GW1 = Q(:,end);
    end
clear Q
end
if int == 1
    reff_max = Y_GW1;
    clear Y_GW1
else
    muW1(:,int) = Y_GW1;
    muW2(:,int) = Y_GW2;
    muW3(:,int) = Y_GW3;
end
%
xdt(int)=log10(dt);  % logarithm of time step wrt base '10'
clear Y_GW1 Y_GW2
%
end
%
%% Plotting of the error
for int = 2:11
errW1(:,int) = log10(norm (muW1(:,int) -  reff_max));
errW2(:,int) = log10(norm (muW2(:,int) -  reff_max));
errW3(:,int) = log10(norm (muW3(:,int) -  reff_max));
%
end
%
ID = [5 6 7 8 9 10 11];
%
YMatrix1 = [errW1(ID);errW2(ID);errW3(ID)];
%
createfigure_fig5(xdt(ID), YMatrix1)
%
% End