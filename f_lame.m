function [o1, o2] = f_lame(p1, p2, fromStr, toStr)
    switch upper(fromStr)
        case "E_NU"
            E = p1; v = p2;
            switch upper(toStr)
                case "LAMBDA_MU"
                    lambda = E*v/(1+v)/(1-2*v);
                    mu = E/2/(1+v);
                    o1 = lambda; o2 = mu;
                otherwise
                    disp("no toStr")
            end
        otherwise
            disp("no fromStr")
    end    
end
