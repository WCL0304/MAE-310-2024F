function [u, v, sigma_xx, sigma_yy, sigma_xy] = manufactured_solution(x, y, L, R, T)
    % 生成制造解
    r = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    u = T * (1 - R^2 / r^2) * cos(theta);
    v = T * (1 - R^2 / r^2) * sin(theta);
    sigma_xx = T * (1 - R^2 / r^2) * (1 + 4 * (R^2 / r^2) - 3 * (R^4 / r^4)) * cos(2 * theta);
    sigma_yy = T * (1 - R^2 / r^2) * (1 - (R^2 / r^2) + 3 * (R^4 / r^4)) * cos(2 * theta);
    sigma_xy = -T * (1 - R^2 / r^2) * (2 * (R^2 / r^2) - 3 * (R^4 / r^4)) * sin(2 * theta);
end