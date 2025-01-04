function [Na_xi, Na_eta] = Quad_grad_tri(aa, xi, eta)
    if aa == 1
        Na_xi = -1;
        Na_eta = -1;
    elseif aa == 2
        Na_xi = 1;
        Na_eta = 0;
    elseif aa == 3
        Na_xi = 0;
        Na_eta = 1;
    else
        error('Error: value of a should be 1, 2, or 3.');
    end
end