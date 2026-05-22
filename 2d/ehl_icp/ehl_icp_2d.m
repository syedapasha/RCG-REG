clear all; close all; 

% shape = "lizard";
shape = "dog";


%% Load source and target shapes
%*******************************
srcFile = strcat('../data/',shape,'-s.mat');
tgtFile = strcat('../data/',shape,'-t.mat');

S = load(srcFile); X1 = S.contour;  % source     
T = load(tgtFile); X2 = T.contour;  % target     

N = min(size(X1,1), size(X2,1));


eta = .5;           % step size
mu = 0.35;       
epsilon = 1e-5;     % stopping threshold
max_iter = 400; 

errs = zeros(2,1);
results = cell(2,1);

for flipFlag = 0:1
    if flipFlag
        X1_try = flipud(X1);
    else
        X1_try = X1;
    end

   
    x_c = mean(X1_try,1)';
    y_c = mean(X2,1)';
    X1_c = X1_try - x_c'; % N x 2


    R = eye(2);
    t = y_c - R * x_c;      
    ER = zeros(max_iter,1);
    J_prev = zeros(2); 


    for k = 1:max_iter
        
        X1_transformed = (R * X1_c') + t; X1_transformed = X1_transformed';
       
        Z = zeros(size(X1_transformed));
        for i = 1:N
            diffs = X2 - X1_transformed(i,:);
            [~, idx] = min(sum(diffs.^2,2));
            Z(i,:) = X2(idx,:);
        end

        ER(k) = sqrt(mean(sum((X1_transformed - Z).^2,2)));
        
%         z_c = mean(Z,1)';
%         Zc = Z - z_c';

        % Step 7: Compute gradient wrt rotation (Eq. 26)
        G = zeros(2);
        for i = 1:N
            xi = X1_c(i,:)';      
%             zi = Zc(i,:)';        
%             G = G + xi * (R*xi + t - z_c - zi)'; 
            zi = Z(i,:)';        
            G = G + xi * (R*xi + t - zi)'; 
        end
        G = (2/N) * G;

        % Step 8: Skew-symmetrize (project onto so(2))
        dR = R' * G;           % Pull back to tangent space
        dR = 0.5 * (dR - dR'); % Skew-symmetric part only

        % Step 9: Hamiltonian update via exponential map with momentum
        J = mu * J_prev + eta * dR;
        R = R * expm(J);
        J_prev = J;

        % Step 10: Update translation (Eq. 24)
        t = y_c - R * x_c;

        % Step 11: Convergence check
        if k > 1 && abs(ER(k-1) - ER(k)) < epsilon
            ER = ER(1:k); % Trim error vector
            break;
        end
    end

    % Store result
    errs(flipFlag+1) = ER(end);
    results{flipFlag+1} = struct('R', R, 't', t, ...
        'Xreg', (R * X1_try' + t)', 'ER', ER, 'k', k, 'X1_try', X1_try);
end

% Choose best (lowest error)
[~, bestIdx] = min(errs);
best = results{bestIdx};

Xreg = best.Xreg;
R = best.R;
t = best.t;
ER = best.ER;
k = best.k;
X1_used = best.X1_try;

if shape == "lizard"
    lizard_ER_ehl_icp = ER;
    save('lizard_ER_ehl_icp.mat', "lizard_ER_ehl_icp");
elseif shape == "dog"
    dog_ER_ehl_icp = ER;
    save('dog_ER_ehl_icp.mat', "dog_ER_ehl_icp");
end    


%% construct plots
%*****************

% plot registered shape
f = figure('visible', 'off');
plot(X2(:,1), X2(:,2), Xreg(:,1), Xreg(:,2), '-.', 'MarkerSize',2);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('EHL-ICP', 'FontSize',15);
legend('target', 'registered', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_ehl_icp.png'));


%% --- Console Output ---
fprintf('\n--- EHL-ICP Hamiltonian Registration (Paper Equations, Auto-Flip, Fast) ---\n');
if bestIdx == 2
    fprintf('Flipped source shape was used for best registration.\n');
else
    fprintf('Original source shape was used for best registration.\n');
end
fprintf('Iterations: %d\n', k);
fprintf('Final RMS Error: %.5f\n', ER(end));
