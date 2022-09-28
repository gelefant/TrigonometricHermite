function Sketch(S,D)
% -------------------------------------------------------------------------
% In this demo we compute the Hermite interpolant for Berrut's barycentric
% trigonometric interpolant 
%
% The demo allows to reconstruct a sketch starting from some extracted
% nodes, in particular it is possibile to increase the value of the integer
% S to use a portion of the initial nodes. 
% 
% INPUT:
% - S: Subdivision of the set of nodes extracted from the sketch, i.e.,
%      #NODES/S
% - D: Different Draws 1 - Dog, 2 - Elephant, 3 - Scrooge McDuck
% 
% 
% -------------------------------------------------------------------------
% Dates
%--------------------------------------------------------------------------
% First version: June 20, 2022;
% Checked: September 28, 2022.
%--------------------------------------------------------------------------
% Authors
%--------------------------------------------------------------------------
% G. Elefante
%--------------------------------------------------------------------------
% Paper
%--------------------------------------------------------------------------
% "A barycentric trigonometric Hermite interpolant via an iterative 
% approach"
% G. Elefante
%--------------------------------------------------------------------------

switch D
    case 1
        T = xlsread('Doggo.xls');
    case 2
        T = xlsread('Elephant.xls');
    case 3
        T = xlsread('ScroogeMcDuck2.xls');
end

N1 = size(T,1);

T1 = [T; T(1,1), T(1,2)];

V = (T1(2:N1+1,:)-T1(1:N1,:))/(2*pi/(N1+1));

T = T(1:S:N1,:);
V = V(1:S:N1,:);
N = size(T,1);

nodes = linspace(0,2*pi,N+1);
nodes(end) = [];

XX = (nodes'-nodes)/2;
XX = XX + eye(size(XX));
id = 0:(N-1);
ID = id-id';

x_eval = linspace(0,2*pi,800);

% X coord
j = 0;
D1_0 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
D1_0 = D1_0 - diag(diag(D1_0));
D1_0 = D1_0 - diag(sum(D1_0,2));
Dx = D1_0*T(:,1);

d_x = @(y,i) 2*sin((y-nodes(i))/2);

b_x = @(y,i) (-1)^(i-1)*cst((y-nodes(i))/2,N)/sum((-1).^id.*cst((y-nodes)/2,N));



for j =1 :length(x_eval)
    Ix(j) = 0;
    for s = 0:1 %Order derivatives 
        for l =1:length(nodes)
            if s == 0
                g = T(:,1);
            elseif s == 1
                g = V(:,1) - Dx;
            end
            Ix(j) = Ix(j) + 1/factorial(s)*d_x(x_eval(j),l)^s*b_x(x_eval(j),l)^(s+1)*g(l);
        end
    end
end

li = ismember(x_eval,nodes);
liN = ismember(nodes,x_eval(li));
Ix(li) = T(liN,1);

% X coord
j = 0;
D1_0 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
D1_0 = D1_0 - diag(diag(D1_0));
D1_0 = D1_0 - diag(sum(D1_0,2));
Dy = D1_0*T(:,2);

d_y = @(y,i) 2*sin((y-nodes(i))/2);

b_y = @(y,i) (-1)^(i-1)*cst((y-nodes(i))/2,N)/sum((-1).^id.*cst((y-nodes)/2,N));

for j =1 :length(x_eval)
    Iy(j) = 0;
    for s = 0:1 %Order derivatives 
        for l =1:length(nodes)
            if s == 0
                g = T(:,2);
            elseif s == 1
                g = V(:,2) - Dy;
            end
            Iy(j) = Iy(j) + 1/factorial(s)*d_y(x_eval(j),l)^s*b_y(x_eval(j),l)^(s+1)*g(l);
        end
    end
end

li = ismember(x_eval,nodes);
liN = ismember(nodes,x_eval(li));
Iy(li) = T(liN,2);

figure(5)
plot(Ix,Iy,'b')
hold on
plot(T(:,1),T(:,2),'rd','MarkerSize',2)

end