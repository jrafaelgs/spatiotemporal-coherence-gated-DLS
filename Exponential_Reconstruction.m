function Exp_Reconst = Exponential_Reconstruction(t,exp_options,Fit_Params)

switch exp_options
    case 1  % Reconstruction of 1 Exponential
        a1 = Fit_Params(1);
        Tau1 = Fit_Params(2);

        expSum_1_1L = (a1*exp(-t*2*pi/Tau1));
        expSum_1L = expSum_1_1L;
        
        Exp_Reconst = [expSum_1_1L expSum_1L];
        
    case 2  % Reconstruction of 2 Exponentials
        a1 = Fit_Params(1);     a2 = Fit_Params(2);
        Tau1 = Fit_Params(3);   Tau2 = Fit_Params(4);

        expSum_1_2L = (a1*exp(-t*2*pi/Tau1));
        expSum_2_2L = (a2*exp(-t*2*pi/Tau2));
        expSum_2L = expSum_1_2L + expSum_2_2L;
        
        Exp_Reconst = [expSum_1_2L expSum_2_2L expSum_2L];
        
    case 3  % Reconstruction of 3 Exponentials
        a1 = Fit_Params(1);     a2 = Fit_Params(2);     a3 = Fit_Params(3);
        Tau1 = Fit_Params(4);   Tau2 = Fit_Params(5);   Tau3 = Fit_Params(6);

        expSum_1_3L = (a1*exp(-t*2*pi/Tau1));
        expSum_2_3L = (a2*exp(-t*2*pi/Tau2));
        expSum_3_3L = (a3*exp(-t*2*pi/Tau3));
        expSum_3L = expSum_1_3L + expSum_2_3L + expSum_3_3L;
        
        Exp_Reconst = [expSum_1_3L expSum_2_3L expSum_3_3L expSum_3L];

    case 4  % Reconstruction of 4 Exponentials
        a1 = Fit_Params(1);     a2 = Fit_Params(2);
        a3 = Fit_Params(3);     a4 = Fit_Params(4);
        Tau1 = Fit_Params(5);   Tau2 = Fit_Params(6);
        Tau3 = Fit_Params(7);   Tau4 = Fit_Params(8);

        expSum_1_4L = (a1*exp(-t*2*pi/Tau1));
        expSum_2_4L = (a2*exp(-t*2*pi/Tau2));
        expSum_3_4L = (a3*exp(-t*2*pi/Tau3));
        expSum_4_4L = (a4*exp(-t*2*pi/Tau4));
        expSum_4L = expSum_1_4L + expSum_2_4L + expSum_3_4L + expSum_4_4L;
        
        Exp_Reconst = [expSum_1_4L expSum_2_4L expSum_3_4L expSum_4_4L expSum_4L];

    case 5  % Reconstruction of 5 Exponentials
        a1 = Fit_Params(1);     a2 = Fit_Params(2);
        a3 = Fit_Params(3);     a4 = Fit_Params(4);
        a5 = Fit_Params(5);
        Tau1 = Fit_Params(6);   Tau2 = Fit_Params(7);
        Tau3 = Fit_Params(8);   Tau4 = Fit_Params(9);
        Tau5 = Fit_Params(10);

        expSum_1_5L = (a1*exp(-t*2*pi/Tau1));
        expSum_2_5L = (a2*exp(-t*2*pi/Tau2));
        expSum_3_5L = (a3*exp(-t*2*pi/Tau3));
        expSum_4_5L = (a4*exp(-t*2*pi/Tau4));
        expSum_5_5L = (a5*exp(-t*2*pi/Tau5));
        expSum_5L = expSum_1_5L + expSum_2_5L + expSum_3_5L + expSum_4_5L + expSum_5_5L;
        
        Exp_Reconst = [expSum_1_5L expSum_2_5L expSum_3_5L expSum_4_5L expSum_5_5L expSum_5L];
        
    otherwise
        
end
