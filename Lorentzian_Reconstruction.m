function Lorent_Reconst = Lorentzian_Reconstruction(xFreqn,lorent_options,Fit_Params)

switch lorent_options
    case 1  % Reconstruction of 1 Lorentzian
        multiLorentz_1L = @(f,x) 2/pi*((x(1)/x(2))./(f.^2 + 1/x(2)^2));
        err_1L = @(f,x) sum((log(multiLorentz_1L(f,x)) - log(PSDn)).^2);

        a1 = Fit_Params(1);
        Tau1 = Fit_Params(2);

        Lorentzian_1_1L = (2/pi)*(a1/Tau1)./(xFreqn.^2 + 1/Tau1^2);
        PSDn_Fit_1L = multiLorentz_1L(xFreqn,[a1 Tau1]);
        
        Lorent_Reconst = [Lorentzian_1_1L PSDn_Fit_1L];
        
    case 2  % Reconstruction of 2 Lorentzians
        multiLorentz_2L = @(f,x) 2/pi*((x(1)/x(3))./(f.^2 + 1/x(3)^2) + (x(2)/x(4))./(f.^2 + 1/x(4)^2));
        err_2L = @(f,x) sum((log(multiLorentz_2L(f,x)) - log(PSDn)).^2);

        a1 = Fit_Params(1);     a2 = Fit_Params(2);
        Tau1 = Fit_Params(3);   Tau2 = Fit_Params(4);

        Lorentzian_1_2L = 2/pi*a1/Tau1./(xFreqn.^2+1/Tau1^2);
        Lorentzian_2_2L = 2/pi*a2/Tau2./(xFreqn.^2+1/Tau2^2);
        PSDn_Fit_2L = multiLorentz_2L(xFreqn,[a1 a2 Tau1 Tau2]);

        Lorent_Reconst = [Lorentzian_1_2L Lorentzian_2_2L PSDn_Fit_2L];
        
    case 3  % Reconstruction of 3 Lorentzians
        multiLorentz_3L = @(f,x) 2/pi*((x(1)/x(4))./(f.^2 + 1/x(4)^2) + (x(2)/x(5))./(f.^2 + 1/x(5)^2) + (x(3)/x(6))./(f.^2 + 1/x(6)^2));
        err_3L = @(f,x) sum((log(multiLorentz_3L(f,x)) - log(PSDn)).^2);

        a1 = Fit_Params(1);     a2 = Fit_Params(2);     a3 = Fit_Params(3);
        Tau1 = Fit_Params(4);   Tau2 = Fit_Params(5);   Tau3 = Fit_Params(6);

        Lorentzian_1_3L = 2/pi*a1/Tau1./(xFreqn.^2+1/Tau1^2);
        Lorentzian_2_3L = 2/pi*a2/Tau2./(xFreqn.^2+1/Tau2^2);
        Lorentzian_3_3L = 2/pi*a3/Tau3./(xFreqn.^2+1/Tau3^2);
        PSDn_Fit_3L = multiLorentz_3L(xFreqn,[a1 a2 a3 Tau1 Tau2 Tau3]);
        
        Lorent_Reconst = [Lorentzian_1_3L Lorentzian_2_3L Lorentzian_3_3L PSDn_Fit_3L];

    case 4  % Reconstruction of 4 Lorentzians
        multiLorentz_4L = @(f,x) 2/pi*((x(1)/x(5))./(f.^2 + 1/x(5)^2) + (x(2)/x(6))./(f.^2 + 1/x(6)^2) + (x(3)/x(7))./(f.^2 + 1/x(7)^2) + (x(4)/x(8))./(f.^2 + 1/x(8)^2));
        err_4L = @(f,x) sum((log(multiLorentz_4L(f,x)) - log(PSDn)).^2);

        a1 = Fit_Params(1);     a2 = Fit_Params(2);     a3 = Fit_Params(3);     a4 = Fit_Params(4);
        Tau1 = Fit_Params(5);   Tau2 = Fit_Params(6);   Tau3 = Fit_Params(7);   Tau4 = Fit_Params(8);

        Lorentzian_1_4L = 2/pi*a1/Tau1./(xFreqn.^2+1/Tau1^2);
        Lorentzian_2_4L = 2/pi*a2/Tau2./(xFreqn.^2+1/Tau2^2);
        Lorentzian_3_4L = 2/pi*a3/Tau3./(xFreqn.^2+1/Tau3^2);
        Lorentzian_4_4L = 2/pi*a4/Tau4./(xFreqn.^2+1/Tau4^2);
        PSDn_Fit_4L = multiLorentz_4L(xFreqn,[a1 a2 a3 a4 Tau1 Tau2 Tau3 Tau4]);
        
        Lorent_Reconst = [Lorentzian_1_4L Lorentzian_2_4L Lorentzian_3_4L Lorentzian_4_4L PSDn_Fit_4L];

    case 5  % Reconstruction of 5 Lorentzians
        multiLorentz_5L = @(f,x) 2/pi*((x(1)/x(6))./(f.^2 + 1/x(6)^2) + (x(2)/x(7))./(f.^2 + 1/x(7)^2) + (x(3)/x(8))./(f.^2 + 1/x(8)^2) + (x(4)/x(9))./(f.^2 + 1/x(9)^2) + (x(5)/x(10))./(f.^2 + 1/x(10)^2));
        err_5L = @(f,x) sum((log(multiLorentz_5L(f,x)) - log(PSDn)).^2);

        a1 = Fit_Params(1);     a2 = Fit_Params(2);     a3 = Fit_Params(3);
        a4 = Fit_Params(4);     a5 = Fit_Params(5);
        Tau1 = Fit_Params(6);   Tau2 = Fit_Params(7);   Tau3 = Fit_Params(8);
        Tau4 = Fit_Params(9);   Tau5 = Fit_Params(10);

        Lorentzian_1_5L = 2/pi*a1/Tau1./(xFreqn.^2+1/Tau1^2);
        Lorentzian_2_5L = 2/pi*a2/Tau2./(xFreqn.^2+1/Tau2^2);
        Lorentzian_3_5L = 2/pi*a3/Tau3./(xFreqn.^2+1/Tau3^2);
        Lorentzian_4_5L = 2/pi*a4/Tau4./(xFreqn.^2+1/Tau4^2);
        Lorentzian_5_5L = 2/pi*a5/Tau5./(xFreqn.^2+1/Tau5^2);
        PSDn_Fit_5L = multiLorentz_5L(xFreqn,[a1 a2 a3 a4 a5 Tau1 Tau2 Tau3 Tau4 Tau5]);
        
        Lorent_Reconst = [Lorentzian_1_5L Lorentzian_2_5L Lorentzian_3_5L Lorentzian_4_5L Lorentzian_5_5L PSDn_Fit_5L];
        
    otherwise
        
end
