function Fit_Params = Lorentzian_Fit(xFreqn,PSDn,lorent_options,x0,random_fit,Ntimes)

min_a = 0.0;        max_a = 1.0;
min_Tau = 0.0;      max_Tau = 1.0;

switch lorent_options
    case 1  % Fit with 1 Lorentzian
        multiLorentz_1L = @(f,x) 2/pi*((x(1)/x(2))./(f.^2 + 1/x(2)^2));
        err_1L = @(f,x) sum((log(multiLorentz_1L(f,x)) - log(PSDn)).^2);

        % Definitions and bounds of fit parameters
        minCharTime = 1/max(xFreqn);
        maxCharTime = 1/min(xFreqn);
        Aeq = vertcat([1 1],zeros(3,2));
        beq = [1 0 0 0]';
        A = [];
        b = [];
        nonlcon = [];
        lb = zeros(1,2);
        ub = ones(1,2);
        options = optimset('Algorithm','sqp');

        for i = 1:Ntimes
            if random_fit == 1
                x0 = [rand(1) rand(1)];  % Random initial guess
                a1_0 = min_a + (max_a - min_a)*x0(1,1);
                Tau1_0 = min_Tau + (max_Tau - min_Tau)*x0(1,2);
            else
                a1_0 = x0(1,1);
                Tau1_0 = x0(1,2);
            end
            x0_vec_1L(i,:) = [a1_0 Tau1_0];    % Initial Guess

            % Perform a multi-Lorentzian fit to the log of the PSD data.
            fitParam_vec_1L(i,:) = fmincon(@(x) err_1L(xFreqn,x), x0_vec_1L(i,:), A, b, Aeq, beq, lb, ub, nonlcon, options);
        end
        
        %  Rounded version of the Fit Parameters
        fitParam_vec_1L(:,1) = round(fitParam_vec_1L(:,1)*(1e4))/1e4;     % Round to 4 decimals
        fitParam_vec_1L(:,2) = round(fitParam_vec_1L(:,2)*(1e6))/1e6;     % Round to 6 decimals

        % Statistics calculated from the ROUNDED fit parameters
        [a1_mode a1_mode_occur]= mode(fitParam_vec_1L(:,1));
        [Tau1_mode Tau1_mode_occur] = mode(fitParam_vec_1L(:,2));

        %  Arbitrarily choosing the mode as the best fit option
        a1 = a1_mode;   Tau1 = Tau1_mode;
        a1_occur = a1_mode_occur;   Tau1_occur = Tau1_mode_occur;
        clc
        fprintf('%1.6g\t %1.6g\t \n%d\t %d\t \n',a1,Tau1,a1_occur,Tau1_occur)

        Fit_Params = [a1;Tau1];
        
    case 2  % Fit with 2 Lorentzians
        multiLorentz_2L = @(f,x) 2/pi*((x(1)/x(3))./(f.^2 + 1/x(3)^2) + (x(2)/x(4))./(f.^2 + 1/x(4)^2));
        err_2L = @(f,x) sum((log(multiLorentz_2L(f,x)) - log(PSDn)).^2);

        % Definitions and bounds of fit parameters
        minCharTime = 1/max(xFreqn);
        maxCharTime = 1/min(xFreqn);
        Aeq = vertcat([1 1 0 0],zeros(3,4));
        beq = [1 0 0 0]';
        A = [];
        b = [];
        nonlcon = [];
        lb = zeros(1,4);
        ub = ones(1,4);
        options = optimset('Algorithm','sqp');

        for i = 1:Ntimes
            if random_fit == 1
                x0 = [rand(1) rand(1) rand(1) rand(1)];  % Random initial guess
                a1_0 = min_a + (max_a - min_a)*x0(1,1);
                a2_0 = min_a + (max_a - min_a)*x0(1,2);
                Tau1_0 = min_Tau + (max_Tau - min_Tau)*x0(1,3);
                Tau2_0 = min_Tau + (max_Tau - min_Tau)*x0(1,4);
            else
                a1_0 = x0(1,1);     a2_0 = x0(1,2);
                Tau1_0 = x0(1,3);   Tau2_0 = min_Tau + (max_Tau - min_Tau)*x0(1,4);
            end
            x0_vec_2L(i,:) = [a1_0 a2_0 Tau1_0 Tau2_0];    % Initial guess

            % Perform a multi-Lorentzian fit to the log of the PSD data.
            fitParam_vec_2L(i,:) = fmincon(@(x) err_2L(xFreqn,x), x0_vec_2L(i,:), A, b, Aeq, beq, lb, ub, nonlcon, options);

            % Forcing the algorithm to give the fitParams in decending order!
            if fitParam_vec_2L(i,4) > fitParam_vec_2L(i,3)
                aux = [fitParam_vec_2L(i,2) fitParam_vec_2L(i,1) fitParam_vec_2L(i,4) fitParam_vec_2L(i,3)];
                fitParam_vec_2L(i,:) = aux;
            end
        end

        %  Rounded version of the Fit Parameters
        fitParam_vec_2L(:,1:2) = round(fitParam_vec_2L(:,1:2)*(1e4))/1e4;     % Rounded to 4 decimals
        fitParam_vec_2L(:,3:4) = round(fitParam_vec_2L(:,3:4)*(1e6))/1e6;     % Rounded to 6 decimals

        % Statistics calculated from the ROUNDED fit parameters
        [a1_mode a1_mode_occur] = mode(fitParam_vec_2L(:,1));
        [a2_mode a2_mode_occur] = mode(fitParam_vec_2L(:,2));
        [Tau1_mode Tau1_mode_occur] = mode(fitParam_vec_2L(:,3));
        [Tau2_mode Tau2_mode_occur] = mode(fitParam_vec_2L(:,4));

        a1 = a1_mode;   a2 = a2_mode;
        a1_occur = a1_mode_occur;   a2_occur = a2_mode_occur;

        Tau1 = Tau1_mode;   Tau2 = Tau2_mode;
        Tau1_occur = Tau1_mode_occur;   Tau2_occur = Tau2_mode_occur;

        clc
        fprintf('%1.6g\t %1.6g\t %1.6g\t %1.6g\n',a1,a2,Tau1,Tau2)
        fprintf('%d\t %d\t %d\t %d\n',a1_occur,a2_occur,Tau1_occur,Tau2_occur)

        Fit_Params = [a1;a2;Tau1;Tau2];

    case 3  % Fit with 3 Lorentzians
        multiLorentz_3L = @(f,x) 2/pi*((x(1)/x(4))./(f.^2 + 1/x(4)^2) + (x(2)/x(5))./(f.^2 + 1/x(5)^2) + (x(3)/x(6))./(f.^2 + 1/x(6)^2));
        err_3L = @(f,x) sum((log(multiLorentz_3L(f,x)) - log(PSDn)).^2);

        % Definitions and bounds of fit parameters
        minCharTime = 1/max(xFreqn);
        maxCharTime = 1/min(xFreqn);
        Aeq = vertcat([1 1 1 0 0 0],zeros(5,6));
        beq = [1 0 0 0 0 0]';
        A = [];
        b = [];
        nonlcon = [];
        lb = zeros(1,6);
        ub = ones(1,6);
        options = optimset('Algorithm','sqp');

        for i = 1:Ntimes
            if random_fit == 1
                x0 = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];  % Random initial guess
                a1_0 = min_a + (max_a - min_a)*x0(1,1);
                a2_0 = min_a + (max_a - min_a)*x0(1,2);
                a3_0 = min_a + (max_a - min_a)*x0(1,3);
                Tau1_0 = min_Tau + (max_Tau - min_Tau)*x0(1,4);
                Tau2_0 = min_Tau + (max_Tau - min_Tau)*x0(1,5);
                Tau3_0 = min_Tau + (max_Tau - min_Tau)*x0(1,6);
            else
                a1_0 = x0(1,1);     a2_0 = x0(1,2);     a3_0 = x0(1,3);
                Tau1_0 = x0(1,4);   Tau2_0 = x0(1,5);   Tau3_0 = x0(1,6);
            end
            x0_vec_3L(i,:) = [a1_0 a2_0 a3_0 Tau1_0 Tau2_0 Tau3_0];

            % Perform a multi-Lorentzian fit to the log of the PSD data.
            [fitParam_vec_3L(i,:) , ~, exitflag] = fmincon(@(x) err_3L(xFreqn,x), x0_vec_3L(i,:), A, b, Aeq, beq, lb, ub, nonlcon, options);
            
%             if exitflag
%                 fitParam_vec_3L(i,:) = ones(1,6);
% %                 return
%             end
            % Forcing the algorithm to give the fitParams in decending order!
            aux = fitParam_vec_3L(i,:);
            if aux(1,6) > aux(1,5)
                aux1 = [aux(1,1) aux(1,3) aux(1,2) aux(1,4) aux(1,6) aux(1,5)];
                aux = aux1;
                if aux(1,5) > aux(1,4)
                    aux1 = [aux(1,2) aux(1,1) aux(1,3) aux(1,5) aux(1,4) aux(1,6)];
                    aux = aux1;
                end
            end
            fitParam_vec_3L(i,:) = aux;
        end

        %  Rounded version of the Fit Parameters
        fitParam_vec_3L(:,1:3) = round(fitParam_vec_3L(:,1:3)*(1e4))/1e4;     % Rounded to 4 decimals
        fitParam_vec_3L(:,4:6) = round(fitParam_vec_3L(:,4:6)*(1e6))/1e6;     % Rounded to 6 decimals

        % Statistics calculated from the ROUNDED fit parameters
        [a1_mode a1_mode_occur] = mode(fitParam_vec_3L(:,1));
        [a2_mode a2_mode_occur] = mode(fitParam_vec_3L(:,2));
        [a3_mode a3_mode_occur] = mode(fitParam_vec_3L(:,3));
        [Tau1_mode Tau1_mode_occur] = mode(fitParam_vec_3L(:,4));
        [Tau2_mode Tau2_mode_occur] = mode(fitParam_vec_3L(:,5));
        [Tau3_mode Tau3_mode_occur] = mode(fitParam_vec_3L(:,6));

        a1 = a1_mode;   a2 = a2_mode;   a3 = a3_mode;
        a1_occur = a1_mode_occur;   a2_occur = a2_mode_occur;   a3_occur = a3_mode_occur;

        Tau1 = Tau1_mode;   Tau2 = Tau2_mode;   Tau3 = Tau3_mode;
        Tau1_occur = Tau1_mode_occur;   Tau2_occur = Tau2_mode_occur;   Tau3_occur = Tau3_mode_occur;

        clc
        fprintf('\n%1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\n\n',a1,a2,a3,Tau1,Tau2,Tau3)
        fprintf('%d\t %d\t %d\t %d\t %d\t %d\n',a1_occur,a2_occur,a3_occur,Tau1_occur,Tau2_occur,Tau3_occur)

        Fit_Params = [a1;a2;a3;Tau1;Tau2;Tau3];

    case 4  % Fit with 4 Lorentzians
        multiLorentz_4L = @(f,x) 2/pi*((x(1)/x(5))./(f.^2 + 1/x(5)^2) + (x(2)/x(6))./(f.^2 + 1/x(6)^2) + (x(3)/x(7))./(f.^2 + 1/x(7)^2) + (x(4)/x(8))./(f.^2 + 1/x(8)^2));
        err_4L = @(f,x) sum((log(multiLorentz_4L(f,x)) - log(PSDn)).^2);

        % Definitions and bounds of fit parameters
        minCharTime = 1/max(xFreqn);
        maxCharTime = 1/min(xFreqn);

        Aeq = vertcat([1 1 1 1 0 0 0 0],zeros(7,8));
        beq = [1 0 0 0 0 0 0 0]';
        A = [];
        b = [];
        nonlcon = [];
        lb = zeros(1,8);
        ub = ones(1,8);
        options = optimset('Algorithm','sqp');

        for i = 1:Ntimes
            if random_fit == 1
                x0 = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];  % Random initial guess
                a1_0 = min_a + (max_a - min_a)*x0(1,1);
                a2_0 = min_a + (max_a - min_a)*x0(1,2);
                a3_0 = min_a + (max_a - min_a)*x0(1,3);
                a4_0 = min_a + (max_a - min_a)*x0(1,4);
                Tau1_0 = min_Tau + (max_Tau - min_Tau)*x0(1,5);
                Tau2_0 = min_Tau + (max_Tau - min_Tau)*x0(1,6);
                Tau3_0 = min_Tau + (max_Tau - min_Tau)*x0(1,7);
                Tau4_0 = min_Tau + (max_Tau - min_Tau)*x0(1,8);
            else
                a1_0 = x0(1,1);     a2_0 = x0(1,2);
                a3_0 = x0(1,3);     a4_0 = x0(1,4);
                Tau1_0 = x0(1,5);   Tau2_0 = x0(1,6);
                Tau3_0 = x0(1,7);   Tau4_0 = x0(1,8);
            end
            % Random guess within boundaries
            x0_vec_4L(i,:) = [a1_0 a2_0 a3_0 a4_0 Tau1_0 Tau2_0 Tau3_0 Tau4_0];

            % Perform a multi-Lorentzian fit to the log of the PSD data.
            fitParam_vec_4L(i,:) = fmincon(@(x) err_4L(xFreqn,x), x0_vec_4L(i,:), A, b, Aeq, beq, lb, ub, nonlcon, options);
            
            % Forcing the algorithm to give the fitParams in decending order!
            aux = fitParam_vec_4L(i,:);
            if aux(1,8) > aux(1,7)
                aux1 = [aux(1,1) aux(1,2) aux(1,4) aux(1,3) aux(1,5) aux(1,6) aux(1,8) aux(1,7)];
                aux = aux1;
                if aux(1,7) > aux(1,6)
                    aux1 = [aux(1,1) aux(1,3) aux(1,2) aux(1,4) aux(1,5) aux(1,7) aux(1,6) aux(1,8)];
                    aux = aux1;
                end
            end
            fitParam_vec_4L(i,:) = aux;
        end

        %  Rounded version of the Fit Parameters
        fitParam_vec_4L(:,1:4) = round(fitParam_vec_4L(:,1:4)*(1e4))/1e4;     % Rounded to 4 decimals
        fitParam_vec_4L(:,5:8) = round(fitParam_vec_4L(:,5:8)*(1e6))/1e6;     % Rounded to 6 decimals

        % Statistics calculated from the ROUNDED fit parameters
        [a1_mode a1_mode_occur] = mode(fitParam_vec_4L(:,1));
        [a2_mode a2_mode_occur] = mode(fitParam_vec_4L(:,2));
        [a3_mode a3_mode_occur] = mode(fitParam_vec_4L(:,3));
        [a4_mode a4_mode_occur] = mode(fitParam_vec_4L(:,4));
        [Tau1_mode Tau1_mode_occur] = mode(fitParam_vec_4L(:,5));
        [Tau2_mode Tau2_mode_occur] = mode(fitParam_vec_4L(:,6));
        [Tau3_mode Tau3_mode_occur] = mode(fitParam_vec_4L(:,7));
        [Tau4_mode Tau4_mode_occur] = mode(fitParam_vec_4L(:,8));

        a1 = a1_mode;   a2 = a2_mode;   a3 = a3_mode;   a4 = a4_mode;
        a1_occur = a1_mode_occur;   a2_occur = a2_mode_occur;
        a3_occur = a3_mode_occur;   a4_occur = a4_mode_occur;

        Tau1 = Tau1_mode;   Tau2 = Tau2_mode;   Tau3 = Tau3_mode;   Tau4 = Tau4_mode;
        Tau1_occur = Tau1_mode_occur;   Tau2_occur = Tau2_mode_occur;
        Tau3_occur = Tau3_mode_occur;   Tau4_occur = Tau4_mode_occur;

        clc
        fprintf('\n%1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\n\n',a1,a2,a3,a4,Tau1,Tau2,Tau3,Tau4)
        fprintf('%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n',a1_occur,a2_occur,a3_occur,a4_occur,Tau1_occur,Tau2_occur,Tau3_occur,Tau4_occur)

        Fit_Params = [a1;a2;a3;a4;Tau1;Tau2;Tau3;Tau4];

    case 5  % Fit with 5 Lorentzians
        multiLorentz_5L = @(f,x) 2/pi*((x(1)/x(6))./(f.^2 + 1/x(6)^2) + (x(2)/x(7))./(f.^2 + 1/x(7)^2) + (x(3)/x(8))./(f.^2 + 1/x(8)^2) + (x(4)/x(9))./(f.^2 + 1/x(9)^2) + (x(5)/x(10))./(f.^2 + 1/x(10)^2));
        err_5L = @(f,x) sum((log(multiLorentz_5L(f,x)) - log(PSDn)).^2);

        % Definitions and bounds of fit parameters
        minCharTime = 1/max(xFreqn);
        maxCharTime = 1/min(xFreqn);

        Aeq = vertcat([1 1 1 1 1 0 0 0 0 0],zeros(9,10));
        beq = [1 0 0 0 0 0 0 0 0 0]';
        A = [];
        b = [];
        nonlcon = [];
        lb = zeros(1,10);
        ub = ones(1,10);
        options = optimset('Algorithm','sqp');

        for i = 1:Ntimes
            if random_fit == 1
                x0 = [rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];  % Random initial guess
                a1_0 = min_a + (max_a - min_a)*x0(1,1);
                a2_0 = min_a + (max_a - min_a)*x0(1,2);
                a3_0 = min_a + (max_a - min_a)*x0(1,3);
                a4_0 = min_a + (max_a - min_a)*x0(1,4);
                a5_0 = min_a + (max_a - min_a)*x0(1,5);
                Tau1_0 = min_Tau + (max_Tau - min_Tau)*x0(1,6);
                Tau2_0 = min_Tau + (max_Tau - min_Tau)*x0(1,7);
                Tau3_0 = min_Tau + (max_Tau - min_Tau)*x0(1,8);
                Tau4_0 = min_Tau + (max_Tau - min_Tau)*x0(1,9);
                Tau5_0 = min_Tau + (max_Tau - min_Tau)*x0(1,10);
            else
                a1_0 = x0(1,1);     a2_0 = x0(1,2);     a3_0 = x0(1,3);
                a4_0 = x0(1,4);     a5_0 = x0(1,5);
                Tau1_0 = x0(1,6);   Tau2_0 = x0(1,7);   Tau3_0 = x0(1,8);
                Tau4_0 = x0(1,9);   Tau5_0 = x0(1,10);
            end
            % Random guess within boundaries
            x0_vec_5L(i,:) = [a1_0 a2_0 a3_0 a4_0 a5_0 Tau1_0 Tau2_0 Tau3_0 Tau4_0 Tau5_0];

            % Perform a multi-Lorentzian fit to the log of the PSD data.
            fitParam_vec_5L(i,:) = fmincon(@(x) err_5L(xFreqn,x), x0_vec_5L(i,:), A, b, Aeq, beq, lb, ub, nonlcon, options);
            
            % Forcing the algorithm to give the fitParams in decending order!
%             aux = fitParam_vec_5L(i,:);
%             if aux(1,8) > aux(1,7)
%                 aux1 = [aux(1,1) aux(1,2) aux(1,4) aux(1,3) aux(1,5) aux(1,6) aux(1,8) aux(1,7)];
%                 aux = aux1;
%                 if aux(1,7) > aux(1,6)
%                     aux1 = [aux(1,1) aux(1,3) aux(1,2) aux(1,4) aux(1,5) aux(1,7) aux(1,6) aux(1,8)];
%                     aux = aux1;
%                 end
%             end
%             fitParam_vec_4L(i,:) = aux;
        end

        %  Rounded version of the Fit Parameters
        fitParam_vec_5L(:,1:5) = round(fitParam_vec_5L(:,1:5)*(1e4))/1e4;     % Rounded to 4 decimals
        fitParam_vec_5L(:,6:10) = round(fitParam_vec_5L(:,6:10)*(1e6))/1e6;     % Rounded to 6 decimals

        % Statistics calculated from the ROUNDED fit parameters
        [a1_mode a1_mode_occur] = mode(fitParam_vec_5L(:,1));
        [a2_mode a2_mode_occur] = mode(fitParam_vec_5L(:,2));
        [a3_mode a3_mode_occur] = mode(fitParam_vec_5L(:,3));
        [a4_mode a4_mode_occur] = mode(fitParam_vec_5L(:,4));
        [a5_mode a5_mode_occur] = mode(fitParam_vec_5L(:,5));

        [Tau1_mode Tau1_mode_occur] = mode(fitParam_vec_5L(:,6));
        [Tau2_mode Tau2_mode_occur] = mode(fitParam_vec_5L(:,7));
        [Tau3_mode Tau3_mode_occur] = mode(fitParam_vec_5L(:,8));
        [Tau4_mode Tau4_mode_occur] = mode(fitParam_vec_5L(:,9));
        [Tau5_mode Tau5_mode_occur] = mode(fitParam_vec_5L(:,10));

        a1 = a1_mode;   a2 = a2_mode;   a3 = a3_mode;
        a4 = a4_mode;   a5 = a5_mode;
        a1_occur = a1_mode_occur;   a2_occur = a2_mode_occur;
        a3_occur = a3_mode_occur;   a4_occur = a4_mode_occur;
        a5_occur = a5_mode_occur;

        Tau1 = Tau1_mode;   Tau2 = Tau2_mode;   Tau3 = Tau3_mode;
        Tau4 = Tau4_mode;   Tau5 = Tau5_mode;
        Tau1_occur = Tau1_mode_occur;   Tau2_occur = Tau2_mode_occur;
        Tau3_occur = Tau3_mode_occur;   Tau4_occur = Tau4_mode_occur;
        Tau5_occur = Tau5_mode_occur;

        clc
        fprintf('\n%1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\t %1.6g\n\n',a1,a2,a3,a4,a5,Tau1,Tau2,Tau3,Tau4,Tau5)
        fprintf('\n%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n',a1_occur,a2_occur,a3_occur,a4_occur,a5_occur,Tau1_occur,Tau2_occur,Tau3_occur,Tau4_occur,Tau5_occur)

        Fit_Params = [a1;a2;a3;a4;a5;Tau1;Tau2;Tau3;Tau4;Tau5];
    
    otherwise
        
end
