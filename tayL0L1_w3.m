function [L0a, L0b, L1a, L1b, L1L1b, L0L1b, L1L0b, L0L0b, L1L1a, L0L1a, L1L0a, L0L0a] = tayL0L1_w3(a, b, y)
%TAYL0L1_W3 Computes first- and second-order compositions of Kolmogorov operators
%
%   [L0a, L0b, L1a, L1b, L1L1b, L0L1b, L1L0b, L0L0b, L1L1a, L0L1a, L1L0a, L0L0a] =
%       tayL0L1_w3(a, b, y)
%
%   Inputs:
%       a - symbolic drift vector field (La x 1)
%       b - symbolic diffusion matrix (La x Lb)
%       y - symbolic state vector (La x 1)
%
%   Outputs:
%       L0a, L0b, L1a, L1b     - First-order Kolmogorov operator terms
%       L1L1b, L0L1b, L1L0b    - Second-order compositions acting on b
%       L0L0b                  - L0 composed twice on b
%       L1L1a, L0L1a, L1L0a    - Second-order compositions acting on a
%       L0L0a                  - L0 composed twice on a

La = size(a, 1);
Lb = size(b, 2);

%% L0a
misc = sym(zeros(size(a)));
for j = 1:La
    misc = misc + a(j) * diff(a, y(j));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc + 0.5 * b(i,k) * b(j,k) * diff(diff(a, y(i)), y(j));
        end
    end
end
L0a = misc;
fprintf("L0a is created.\n");

%% L0b
misc = sym(zeros(size(b)));
for i = 1:La
    misc = misc + a(i) * diff(b, y(i));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc + 0.5 * b(i,k) * b(j,k) * diff(diff(b, y(i)), y(j));
        end
    end
end
L0b = misc;
fprintf("L0b is created.\n");

%% L1a
misc = sym(zeros(size(b)));
for i = 1:Lb
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i) * diff(a(k), y(j));
        end
    end
end
L1a = misc;
fprintf("L1a is created.\n");

%% L1b
misc = sym(zeros(size(b)));
for i = 1:Lb
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i) * diff(b(k,i), y(j));
        end
    end
end
L1b = misc;
fprintf("L1b is created.\n");

%% L1L1b
misc = sym(zeros(size(b)));
for i = 1:Lb
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i) * diff(L1b(k,i), y(j));
        end
    end
end
L1L1b = misc;
fprintf("L1L1b is created.\n");

%% L1L0b
misc = sym(zeros(size(b)));
for i = 1:Lb
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i) * diff(L0b(k,i), y(j));
        end
    end
end
L1L0b = misc;
fprintf("L1L0b is created.\n");

%% L0L0b
misc = sym(zeros(size(b)));
for i = 1:La
    misc = misc + a(i) * diff(L0b, y(i));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc + 0.5 * b(i,k) * b(j,k) * diff(diff(L0b, y(i)), y(j));
        end
    end
end
L0L0b = misc;
fprintf("L0L0b is created.\n");

%% L0L1b
misc = sym(zeros(size(b)));
for i = 1:La
    misc = misc + a(i) * diff(L1b, y(i));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc + 0.5 * b(i,k) * b(j,k) * diff(diff(L1b, y(i)), y(j));
        end
    end
end
L0L1b = misc;
fprintf("L0L1b is created.\n");

%% L1L1a
misc = sym(zeros(size(b)));
for i = 1:Lb
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i) * diff(L1a(k,i), y(j));
        end
    end
end
L1L1a = misc;
fprintf("L1L1a is created.\n");

%% L1L0a
misc = sym(zeros(size(b)));
for i = 1:Lb
    for k = 1:La
        for j = 1:La
            misc(k,i) = misc(k,i) + b(j,i) * diff(L0a(k,i), y(j));
        end
    end
end
L1L0a = misc;
fprintf("L1L0a is created.\n");

%% L0L0a
misc = sym(zeros(size(b)));
for i = 1:La
    misc = misc + a(i) * diff(L0a, y(i));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc + 0.5 * b(i,k) * b(j,k) * diff(diff(L0a, y(i)), y(j));
        end
    end
end
L0L0a = misc;
fprintf("L0L0a is created.\n");

%% L0L1a
misc = sym(zeros(size(b)));
for i = 1:La
    misc = misc + a(i) * diff(L1a, y(i));
end
for k = 1:Lb
    for i = 1:La
        for j = 1:La
            misc = misc + 0.5 * b(i,k) * b(j,k) * diff(diff(L1a, y(i)), y(j));
        end
    end
end
L0L1a = misc;
fprintf("L0L1a is created.\n");

end
