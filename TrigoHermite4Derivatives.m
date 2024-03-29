function TrigoHermite4Derivatives(part)
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
Test_Fun = 1;

switch Test_Fun
    case 1
        %Test function 1
        f = @(x) exp(2*sin(x)+cos(x));
        Df = @(x) exp(2*sin(x)+cos(x)).*(2*cos(x)-sin(x));
        D2f = @(x) exp(cos(x) + 2*sin(x)).*(4*cos(x).^2 + (-2 + sin(x)).*sin(x) - cos(x).*(1 + 4*sin(x)));
        D3f = @(x) exp(cos(x) + 2*sin(x)).*(-sin(x).^3 + 6*sin(x).^2 + sin(x) + 8*cos(x).^3 - 6*(2*sin(x) + 1).*cos(x).^2 + (6*sin(x).^2 - 9*sin(x) - 2).*cos(x));
        D4f = @(x) exp(cos(x) + 2*sin(x)).*(sin(x).*(sin(x).^3 - 12*sin(x).^2 + 8*sin(x) + 2) + 16*cos(x).^4 - 8*(4*sin(x) + 3).*cos(x).^3 + (24*sin(x).^2 - 24*sin(x) - 13).*cos(x).^2 + (-8*sin(x).^3 + 42*sin(x).^2 + 28*sin(x) + 1).*cos(x));
    case 2
        %Test function 2
        f = @(x) cos(3*x)./cosh(sin(x));
        Df = @(x) -sech(sin(x)).*(3*sin(3*x) + cos(x).*cos(3*x).*tanh(sin(x)));
        D2f = @(x) sech(sin(x)).*(-9*cos(3*x) + cos(x).^2.*cos(3*x).*tanh(sin(x)).^2 - cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + tanh(sin(x)).*(sin(x).*cos(3*x) + 6*sin(3*x).*cos(x)));
        D3f = @(x) sech(sin(x)).*(27*sin(3*x) - 3*sin(3*x).*cos(x).^2.*tanh(sin(x)).^2 + ...
            3*sin(3*x).*cos(x).^2.*sech(sin(x)).^2 - 2*sin(x).*cos(3*x).*cos(x).*tanh(sin(x)).^2 +...
            tanh(sin(x)).*(19*cos(x).*cos(3*x) - 9*sin(x).*sin(3*x)) + 2*sin(x).*cos(3*x).*cos(x).*sech(sin(x)).^2 + ...
            cos(x).*(sin(x).*cos(3*x) + 6*sin(3*x).*cos(x)).*sech(sin(x)).^2 + ...
            4*cos(3*x).*cos(x).^3.*tanh(sin(x)).*sech(sin(x))).^2 - cos(x).*tanh(sin(x)).*sech(sin(x)).*(-9*cos(3*x) + ...
            cos(x).^2.*cos(3*x).*tanh(sin(x)).^2 - cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + ...
            tanh(sin(x)).*(sin(x).*cos(3*x) + 6*sin(3*x).*cos(x)));
        D4f = @(x) -cos(x).^2.*sech(sin(x)).^3.*(-9*cos(3*x) - cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + ...
            (cos(3*x).*sin(x) + 6*cos(x).*sin(3*x)).*tanh(sin(x)) + cos(x).^2.*cos(3*x).*tanh(sin(x)).^2) + ...
            sech(sin(x)).*sin(x).*tanh(sin(x)).*(-9*cos(3*x) - cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + ...
            (cos(3*x).*sin(x) + 6*cos(x).*sin(3*x)).*tanh(sin(x)) + cos(x).^2.*cos(3*x).*tanh(sin(x)).^2) + ...
            cos(x).^2.*sech(sin(x)).*tanh(sin(x)).^2.*(-9*cos(3*x) - cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + ...
            (cos(3*x).*sin(x) + 6*cos(x).*sin(3*x)).*tanh(sin(x)) + cos(x).^2.*cos(3*x).*tanh(sin(x)).^2) - 2*cos(x).*sech(sin(x)).*tanh(sin(x)).*(2*cos(x).*cos(3*x).*sech(sin(x)).^2.*sin(x) + ...
            27*sin(3*x) + 3*cos(x).^2.*sech(sin(x)).^2.*sin(3*x) + cos(x).*sech(sin(x)).^2.*(cos(3*x).*sin(x) + ...
            6*cos(x).*sin(3*x)) + 4*cos(x).^3.*cos(3*x).*sech(sin(x)).^2.*tanh(sin(x)) + (19*cos(x).*cos(3*x) - 9*sin(x).*sin(3*x)).*tanh(sin(x)) - 2*cos(x).*cos(3*x).*sin(x).*tanh(sin(x)).^2 - 3*cos(x).^2.*sin(3*x).*tanh(sin(x)).^2) + ...
            sech(sin(x)).*(81*cos(3*x) + 11*cos(x).^2.*cos(3*x).*sech(sin(x)).^2 + 4*cos(x).^4.*cos(3*x).*sech(sin(x)).^4 - 2*cos(3*x).*sech(sin(x)).^2.*sin(x).^2 - 12*cos(x).*sech(sin(x)).^2.*sin(x).*sin(3*x) - sech(sin(x)).^2.*sin(x).*(cos(3*x).*sin(x) + 6*cos(x).*sin(3*x)) + ...
            2*cos(x).*sech(sin(x)).^2.*(19*cos(x).*cos(3*x) - 9*sin(x).*sin(3*x)) - 20*cos(x).^2.*cos(3*x).*sech(sin(x)).^2.*sin(x).*tanh(sin(x)) - 24*cos(x).^3.*sech(sin(x)).^2.*sin(3*x).*tanh(sin(x)) + ...
            (-46*cos(3*x).*sin(x) - 66*cos(x).*sin(3*x)).*tanh(sin(x)) - 2*cos(x).^2.*sech(sin(x)).^2.*(cos(3*x).*sin(x) + ...
            6*cos(x).*sin(3*x)).*tanh(sin(x)) - 11*cos(x).^2.*cos(3*x).*tanh(sin(x)).^2 - 8*cos(x).^4.*cos(3*x).*sech(sin(x)).^2.*tanh(sin(x)).^2 + ...
            2*cos(3*x).*sin(x).^2.*tanh(sin(x)).^2 + 12*cos(x).*sin(x).*sin(3*x).*tanh(sin(x)).^2);
    case 3
        f = @(x) cos(3*x) + log(cos(x)+1.5);
        Df = @(x) -sin(x)./(1.5 + cos(x)) - 3*sin(3*x);
        D2f = @(x) -cos(x)./(1.5 + cos(x)) - (sin(x).^2)./(1.5 + cos(x)).^2 - 9*cos(3*x);
        D3f = @(x) -(2*sin(x).^3)./(1.5 + cos(x)).^3 + sin(x)./(1.5 + cos(x)) - (3*sin(x).*cos(x))./(1.5 + cos(x)).^2 + 27*sin(3*x);        
        D4f = @(x) -(3*cos(x).^2)./(1.5 + cos(x)).^2 + cos(x)./(1.5 + cos(x)) - (6*sin(x).^4)./(1.5 + cos(x)).^4 + (4*sin(x).^2)./(1.5 + cos(x)).^2 - (12*sin(x).^2.*cos(x))./(1.5 + cos(x)).^3 + 81*cos(3*x);
    case 4
        f = @(x) tanh(50*cos(x+pi/3));
        Df = @(x)  -50*cos(pi/6 - x).*sech(50*sin(pi/6 - x)).^2;
        D2f = @(x) -50*sech(50*sin(pi/6 - x)).^2.*(sin(pi/6 - x) + 100*cos(pi/6 - x).^2.*tanh(50*sin(pi/6 - x)));
        D3f = @(x) 50*cos(pi/6 - x).*sech(50*sin(pi/6 - x)).^2.*(-300*sin(pi/6 - x).*tanh(50*sin(pi/6 - x)) - 10000*cos(pi/6 - x).^2.*tanh(50*sin(pi/6 - x)).^2 +...
            5000*cos(pi/6 - x).^2.*sech(50*sin(pi/6 - x)).^2 + 1);
        D4f = @(x) 50*cos(pi/6 - x).*sech(50*sin(pi/6 - x)).^2.*(-20000*sin(pi/6 - x).*cos(pi/6 - x).*tanh(50*sin(pi/6 - x)).^2 + ...
            300*cos(pi/6 - x).*tanh(50*sin(pi/6 - x)) + 25000*sin(pi/6 - x).*cos(pi/6 - x).*sech(50*sin(pi/6 - x)).^2 + ...
            1500000*cos(pi/6 - x).^3.*tanh(50*sin(pi/6 - x)).*sech(50*sin(pi/6 - x))).^2 + 50*sin(pi/6 - x).*sech(50*sin(pi/6 - x)).^2.*(-300*sin(pi/6 - x).*tanh(50*sin(pi/6 - x)) - 10000*cos(pi/6 - x).^2.*tanh(50*sin(pi/6 - x)).^2 +...
            5000*cos(pi/6 - x).^2.*sech(50*sin(pi/6 - x)).^2 + 1) + 5000*cos(pi/6 - x).^2.*tanh(50*sin(pi/6 - x)).*sech(50*sin(pi/6 - x)).^2.*(-300*sin(pi/6 - x).*tanh(50*sin(pi/6 - x)) - 10000*cos(pi/6 - x).^2.*tanh(50*sin(pi/6 - x)).^2 + ...
            5000*cos(pi/6 - x).^2.*sech(50*sin(pi/6 - x)).^2 + 1);
    case 5
        f = @(x) tanh(sin(x)).*exp(-cos(x));
        Df = @(x) exp(-cos(x)).*(cos(x).*sech(sin(x)).^2 + sin(x).*tanh(sin(x)));
        D2f = @(x) exp(-cos(x)).*((cos(x) + sin(x).^2).*tanh(sin(x)) + sech(sin(x)).^2.*(-sin(x) + ...
            sin(2*x) - 2*cos(x).^2.*tanh(sin(x))));
        D3f = @(x) exp(-cos(x)).*(-2*cos(x).^3.*sech(sin(x)).^4 + sin(x).*(-1 + 3*cos(x) + sin(x).^2).*tanh(sin(x)) + ...
            sech(sin(x)).^2.*(-(cos(x)/4) + 3*cos(2*x) - 3/4*cos(3*x) - 6*(-1 + cos(x)).*cos(x).*sin(x).*tanh(sin(x)) + ...
            4*cos(x).^3.*tanh(sin(x)).^2));
        D4f = @(x) 1/8*exp(-cos(x)).*((-1 + 4*cos(x) + 24*cos(2*x) - 12*cos(3*x) + cos(4*x)).*tanh(sin(x)) + ...
            32*sech(sin(x)).^4.*((3 - 2*cos(x)).*cos(x).^2.*sin(x) + 4*cos(x).^4.*tanh(sin(x))) - 4*sech(sin(x)).^2.*(sin(x) +...
            12*sin(2*x) - 9*sin(3*x) + sin(4*x) + (1 + 6*cos(x) - 14*cos(2*x) + 18*cos(3*x) - 3*cos(4*x)).*tanh(sin(x)) - 16*cos(x).^2.*(-3*sin(x) +...
            sin(2*x)).*tanh(sin(x)).^2 + 16*cos(x).^4.*tanh(sin(x)).^3));
end    


switch part
    case 1
        N = 60;

        nodes = linspace(0,2*pi,N+1);
        nodes(end) = [];

        if shift_nodes
            alpha1 = 0.5;
            alpha2 = 0.5;
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
        F = f(nodes)'; DF = Df(nodes)'; DDF = D2f(nodes)'; DDDF = D3f(nodes)';
        DDDDF = D4f(nodes)';

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

        j = 2;
        D3_2 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
        D3_2 = D3_2 - diag(diag(D3_2));
        D3_2 = D3_2 - diag(sum(D3_2,2));
        DDDR = D3_2*(DDF - DDR);
        DD2R = D2_1^2*(DF-DR);
        D3R = D1_0^3*F;
        
        j = 3;
        D4_3 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
        D4_3 = D4_3 - diag(diag(D4_3));
        D4_3 = D4_3 - diag(sum(D4_3,2));
        DDDDR = D4_3*(DDDF - DDDR);
        DDD2R = D3_2^2*(DDF-DDR);
        DD3R = D2_1^3*(DF-DR);
        D4R = D1_0^4*F;
        
        d = @(y,i) 2*sin((y-nodes(i))/2);

        b = @(y,i) (-1)^(i-1)*cst((y-nodes(i))/2,N)/sum((-1).^id.*cst((y-nodes)/2,N));

        x_eval = linspace(0,2*pi,1500);


        for j =1 :length(x_eval)
            Int(j) = 0;
            for s = 0:4 %Order derivatives 
                for l =1:length(nodes)
                    if s == 0
                        g = F;
                    elseif s == 1
                        g = DF - DR;
                    elseif s == 2
                        g = DDF - D2R - DDR;
                    elseif s == 3
                        g = DDDF - DDDR - DD2R - D3R;
                    else
                        g = DDDDF - DDDDR - DDD2R - DD3R - D4R;
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
        N_value = [5,10,20,40,80,160,320,640,1280];

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
            F = f(nodes)'; DF = Df(nodes)'; DDF = D2f(nodes)'; DDDF = D3f(nodes)';
            DDDDF = D4f(nodes)';
            
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
            
            j = 2;
            D3_2 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
            D3_2 = D3_2 - diag(diag(D3_2));
            D3_2 = D3_2 - diag(sum(D3_2,2));
            DDDR = D3_2*(DDF - DDR);
            DD2R = D2_1^2*(DF-DR);
            D3R = D1_0^3*F;
            
            j = 3;
            D4_3 = (j+1)*(-1).^((j+1)*ID)./2.*cst(XX,N);
            D4_3 = D4_3 - diag(diag(D4_3));
            D4_3 = D4_3 - diag(sum(D4_3,2));
            DDDDR = D4_3*(DDDF - DDDR);
            DDD2R = D3_2^2*(DDF-DDR);
            DD3R = D2_1^3*(DF-DR);
            D4R = D1_0^4*F;
            
            d = @(y,i) 2*sin((y-nodes(i))/2);

            b = @(y,i) (-1)^(i-1)*cst((y-nodes(i))/2,N)/sum((-1).^id.*cst((y-nodes)/2,N));

            x_eval = linspace(0,2*pi,3500);


            for j =1 :length(x_eval)
                Int(j) = 0;
                for s = 0:4 %Order derivatives 
                    for l =1:length(nodes)
                        if s == 0
                            g = F;
                        elseif s == 1
                            g = DF - DR;
                        elseif s == 2
                            g = DDF - D2R - DDR;
                        elseif s == 3
                            g = DDDF - DDDR - DD2R - D3R;
                        else
                            g = DDDDF - DDDDR - DDD2R - DD3R - D4R;
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

        approx_exp = -log2(err(2:end)./err(1:end-1));

        figure(2)
        semilogy(N_value,err)
        title('Absolute Error with N = 5,10,20,40,...,1280')
        xlabel('N')
        ylabel('Absolute Error')

        fprintf('\n ------------------------------\n ')
        for u=1:(length(N_value)-1)
            fprintf(' Approximate ord. conv. N=%d: %1.2f \n ',N_value(u),approx_exp(u))
        end
        fprintf(' ------------------------------\n ')

end



end