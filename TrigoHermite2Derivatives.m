function TrigoHermite2Derivatives(part)
% -------------------------------------------------------------------------
% In this demo we compute the Hermite interpolant for Berrut's barycentric
% trigonometric interpolant 
%
% The demo is split in two part: the first allow to compute for a fixed N
% the interpolant for a test function - 4 available, the second allow to
% compute the error for a test function once the degree is doubled
%
% The results are presented with some statistics and plots
% 
% -------------------------------------------------------------------------
% Dates
%--------------------------------------------------------------------------
% First version: June 01, 2022;
% Checked: January 27, 2023.
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
shift_nodes = 1;
Test_Fun = 4;

switch Test_Fun
    case 1
        %Test function 1
        f = @(x) exp(2*sin(x)+cos(x));
        Df = @(x) exp(2*sin(x)+cos(x)).*(2*cos(x)-sin(x));
        D2f = @(x) exp(cos(x) + 2*sin(x)).*(4*cos(x).^2 + (-2 + sin(x)).*sin(x) - cos(x).*(1 + 4*sin(x)));
    case 2
        %Test function 2
        f = @(x) cos(3*x)./cosh(sin(x));
        Df = @(x) -sech(sin(x)).*(3*sin(3*x) + cos(x).*cos(3*x).*tanh(sin(x)));
        D2f = @(x) sech(sin(x)).*(-9*cos(3*x) + cos(x).^2.*cos(3*x).*tanh(sin(x)).^2 - cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + tanh(sin(x)).*(sin(x).*cos(3*x) + 6*sin(3*x).*cos(x)));
    case 3
        f = @(x) cos(3*x) + log(cos(x)+1.5);
        Df = @(x) -sin(x)./(1.5 + cos(x)) - 3*sin(3*x);
        D2f = @(x) -cos(x)./(1.5 + cos(x)) - (sin(x).^2)./(1.5 + cos(x)).^2 - 9*cos(3*x);
    case 4
        f = @(x) tanh(50*cos(x+pi/3));
        Df = @(x)  -50*cos(pi/6 - x).*sech(50*sin(pi/6 - x)).^2;
        D2f = @(x) -50*sech(50*sin(pi/6 - x)).^2.*(sin(pi/6 - x) + 100*cos(pi/6 - x).^2.*tanh(50*sin(pi/6 - x)));
    case 5
        f = @(x) tanh(sin(x)).*exp(-cos(x));
        Df = @(x) exp(-cos(x)).*(cos(x).*sech(sin(x)).^2 + sin(x).*tanh(sin(x)));
        D2f = @(x) exp(-cos(x)).*((cos(x) + sin(x).^2).*tanh(sin(x)) + sech(sin(x)).^2.*(-sin(x) + ...
            sin(2*x) - 2*cos(x).^2.*tanh(sin(x))));
end       


switch part
    case 1
        N = 120;


        nodes = linspace(0,2*pi,N+1);
        nodes(end) = [];

        if shift_nodes
            alpha1 = 0.75;
            alpha2 = 0.75;
            TTT1=pi/6;
            TTT2=7/6*pi;
            det = @(a,b,x) abs(exp(1i*x)-exp(-1i*x)*a*b).^2-imag(exp(-1i*x)*(a+b)).^2;
            k1 = @(a,b,x) (1i*imag(exp(-1i*x)*(a+b))+(det(a,b,x)).^(1/2))./...
            (exp(-1i*x)-exp(1i*x)*conj(a)*conj(b));
            y1 = mod(real(-1i*log(k1(alpha1*exp(1i*TTT1),alpha2*exp(1i*TTT2),nodes))),2*pi);
            nodes = sort(y1);
        end

        XX = (nodes'-nodes)/2;
        XX = XX + eye(size(XX));
        id = 0:(N-1);
        ID = id-id';

        %
        F = f(nodes)'; DF = Df(nodes)'; DDF = D2f(nodes)';

        j = 0;
        D1_0 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
        D1_0 = D1_0 - diag(diag(D1_0));
        D1_0 = D1_0 - diag(sum(D1_0,2));
        DR = D1_0*F;

        j = 1;
        D2_1 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
        D2_1 = D2_1 - diag(diag(D2_1));
        D2_1 = D2_1 - diag(sum(D2_1,2));
        DDR = D2_1*(DF-DR);
        D2R = D1_0^2*F;

        d = @(y,i) 2*sin((y-nodes(i))/2);

        b = @(y,i) (-1)^(i-1)*cst((y-nodes(i))/2,N)/sum((-1).^id.*cst((y-nodes)/2,N));

        x_eval = linspace(0,2*pi,1500);


        for j =1 :length(x_eval)
            Int(j) = 0;
            for s = 0:2 %Order derivatives 
                for l =1:length(nodes)
                    if s == 0
                        g = F;
                    elseif s == 1
                        g = DF - DR;
                    else
                        g = DDF - D2R - DDR;
                    end
                    Int(j) = Int(j) + 1/factorial(s)*d(x_eval(j),l)^s*b(x_eval(j),l)^(s+1)*g(l);
                end
            end
        end

        li = ismember(x_eval,nodes);
        liN = ismember(nodes,x_eval(li));
        Int(li) = f(nodes(liN));

        MinG = min([f(x_eval),Int]); 
        MaxG = max([f(x_eval),Int]);

        figure(1)
        plot(x_eval,f(x_eval),'--r');
        hold on;
        plot(x_eval,Int,'b');
        plot(nodes,f(nodes),'.g','MarkerSize', 15)
        axis([0 2*pi MinG-1 MaxG+1])
        % axis square

        err_abs = max(abs(Int-f(x_eval)));

        err_rel = max(abs(Int-f(x_eval))./abs(f(x_eval)));

        fprintf('\n\n Number of nodes: %i \n Relative Error : %2.4e \n Absolute Error : %2.4e \n\n', N, err_rel, err_abs)

    case 2
        N_value = [5,10,20,40,80,160];

        for h = 1:length(N_value)
            N = N_value(h);
            nodes = linspace(0,2*pi,N+1);
            nodes(end) = [];
            
            if shift_nodes
                alpha1 = 0.85;
                alpha2 = 0.85;
                TTT1=pi/6;
                TTT2=7/6*pi;
                det = @(a,b,x) abs(exp(1i*x)-exp(-1i*x)*a*b).^2-imag(exp(-1i*x)*(a+b)).^2;
                k1 = @(a,b,x) (1i*imag(exp(-1i*x)*(a+b))+(det(a,b,x)).^(1/2))./...
                (exp(-1i*x)-exp(1i*x)*conj(a)*conj(b));
                y1 = mod(real(-1i*log(k1(alpha1*exp(1i*TTT1),alpha2*exp(1i*TTT2),nodes))),2*pi);
                nodes = sort(y1);
            end

            XX = (nodes'-nodes)/2;
            XX = XX + eye(size(XX));
            id = 0:(N-1);
            ID = id-id';

            %
            F = f(nodes)'; DF = Df(nodes)'; DDF = D2f(nodes)';

            j = 0;
            D1_0 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
            D1_0 = D1_0 - diag(diag(D1_0));
            D1_0 = D1_0 - diag(sum(D1_0,2));
            DR = D1_0*F;

            j = 1;
            D2_1 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
            D2_1 = D2_1 - diag(diag(D2_1));
            D2_1 = D2_1 - diag(sum(D2_1,2));
            DDR = D2_1*(DF-DR);
            D2R = D1_0^2*F;
            
            d = @(y,i) 2*sin((y-nodes(i))/2);

            b = @(y,i) (-1)^(i-1)*cst((y-nodes(i))/2,N)/sum((-1).^id.*cst((y-nodes)/2,N));

            x_eval = linspace(0,2*pi,1500);


            for j =1 :length(x_eval)
                Int(j) = 0;
                for s = 0:2 %Order derivatives 
                    for l =1:length(nodes)
                        if s == 0
                            g = F;
                        elseif s == 1
                            g = DF - DR;
                        else
                            g = DDF - D2R - DDR;
                        end
                        Int(j) = Int(j) + 1/factorial(s)*d(x_eval(j),l)^s*b(x_eval(j),l)^(s+1)*g(l);
                    end
                end
            end

            li = ismember(x_eval,nodes);
            liN = ismember(nodes,x_eval(li));
            Int(li) = f(nodes(liN));

            err(h) = max(abs(f(x_eval)-Int));

        end

        figure(2)
        semilogy(N_value,err)
        title('Absolute Error with N = 5,10,20,40')
        xlabel('N')
        ylabel('Absolute Error')
end

end