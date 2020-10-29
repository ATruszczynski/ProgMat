% problem 1
% clc
% clear
% format
% format compact
% 
% f = [2;-3];
% A = [1, -2; 6, 5; 3, - 1; 1, 3];
% b = [2; 15; 5; 6];
% LB = [0; -inf];
% UB = [inf; 0];
% 
% options = optimset(@linprog);
% %options = optimset(options, 'Display', 'iter', 'Algorithm', 'interior-point');
% options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');
% 
% [x,fval,exitflag,output,lambda] = linprog(-f, A, b, [], [], LB, UB, [], options)

% wykres
% hold on
% grid on
% x1=-2:0.1:2;
% x2=x1;
% [X1, X2] = meshgrid(x1, x2);
% F = f(1) .* X1 + f(2) .* X2;
% [C,h]=contour(X1,X2,F, 'r-');
% clabel(C,h);
% 
% for i=1:4
%    G = A(i,1) .* X1 + A(i,2) .* X2 - b(i);
%    contour(X1, X2, G, 0:-0.1:-0.5, 'g-');
%    contour(X1, X2, G, [0,0], 'b-');
%    gtext(sprintf('g%d', i));
%    
% end

%cośtam
% x_opt = [A(1,:); A(3,:)]\[b(1); b(3)];
% disp(x_opt);
% f_opt=f'*x_opt;
% disp(f_opt);

% problem 2

clc
clear
format
format compact

f = [5;2;3;5];
A = [-1, 1, -7, -3; 2, 3, 1, -4];
b = [-4; 5];
Aeq = [1, 2, 2, 1];
beq = [9];
LB = [0; 0; 0; 0];
UB = [];

options = optimset(@linprog);
%options = optimset(options, 'Display', 'iter', 'Algorithm', 'interior-point');
options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');

[x,fval,exitflag,output,lambda] = linprog(f, A, b, Aeq, beq, LB, UB, [], options)

% wykres
% hold on
% grid on
% x1=-2:0.1:2;
% x2=x1;
% [X1, X2] = meshgrid(x1, x2);
% F = f(1) .* X1 + f(2) .* X2;
% [C,h]=contour(X1,X2,F, 'r-');
% clabel(C,h);
% 
% for i=1:4
%    G = A(i,1) .* X1 + A(i,2) .* X2 - b(i);
%    contour(X1, X2, G, 0:-0.1:-0.5, 'g-');
%    contour(X1, X2, G, [0,0], 'b-');
%    gtext(sprintf('g%d', i));
%    
% end

%cośtam
% x_opt = [A(1,:); A(3,:)]\[b(1); b(3)];
% disp(x_opt);
% f_opt=f'*x_opt;
% disp(f_opt);

% problem 3

% clc
% clear
% format
% format compact
% 
% f = [10;0;20];
% A = [-2, -1, 0; -1, -4, -6];
% b = [-21; -14];
% Aeq = [];
% beq = [];
% LB = [0; 0; 0];
% UB = [];
% 
% options = optimset(@linprog);
% %options = optimset(options, 'Display', 'iter', 'Algorithm', 'interior-point');
% options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');
% 
% [x,fval,exitflag,output,lambda] = linprog(f, A, b, Aeq, beq, LB, UB, [], options)


% problem 4

clc
clear
format
format compact

f = [1; 1; 1; 1; 1; 1; 1];
A = [1, 0, 0, 1, 1, 1, 1; 1, 1, 0, 0, 1, 1, 1; 1, 1, 1, 0, 0, 1, 1; 1, 1, 1, 1, 0, 0, 1; 1, 1, 1, 1, 1, 0, 0; 0, 1, 1, 1, 1, 1, 0; 0, 0, 1, 1, 1, 1, 1];
A = -1 .* A;
b = [15, 10, 15, 20, 15, 16, 10];
b = -1 .* b;
Aeq = [];
beq = [];
LB = [0; 0; 0; 0; 0; 0; 0];
UB = [];

options = optimset(@linprog);
%options = optimset(options, 'Display', 'iter', 'Algorithm', 'interior-point');
options = optimset(options, 'Display', 'iter', 'Algorithm', 'dual-simplex');

[x,fval,exitflag,output,lambda] = linprog(f, A, b, Aeq, beq, LB, UB, [], options)




















