function [err, grad_R, grad_t] = registration_error_se3(R, t, X1_c, Z)
% Registration error and gradients for SE(3)

% N = size(X1_c,1);
X1_transformed = (R * X1_c') + t; X1_transformed = X1_transformed';
err = sqrt(mean(sum((X1_transformed - Z).^2,2)));

% z_c = mean(Z,1)';
% Zc = Z - z_c';
% 
% G_rot = zeros(3,1); G_t = zeros(3,1);
% for i = 1:N
%     xi = X1_c(i,:)';
%     zi = Zc(i,:)';
%     diff = R*xi + t - z_c - zi;
%     G_rot = G_rot + cross(xi, diff);
%     G_t = G_t + diff;
% end
% grad_R = (2/N) * G_rot;
% grad_t = (2/N) * G_t;


%% for fast computation the for loop is replaced
%
diff = R*X1_c' + t - Z';
grad_R = 2 * mean(cross(X1_c',diff),2);
grad_t = 2 * mean(diff,2);
end