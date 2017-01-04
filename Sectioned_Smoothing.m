function  [psd_Smooth error_Smooth] = Sectioned_Smoothing(xFreq,psd_raw,Breaking_freq,SM_Method_a,Span_a,SM_Method_b,Span_b)

% Breaking the PSD into 2 different regimes for independent smooting
psd_raw_a = psd_raw(xFreq<=Breaking_freq);
psd_raw_b = psd_raw(xFreq>Breaking_freq);

% Smoothing each section of the raw PSD
psd_Smooth_a = smooth(psd_raw_a,Span_a,SM_Method_a);
error_Smooth_a = (abs(psd_raw_a - psd_Smooth_a)).^2;

psd_Smooth_b = smooth(psd_raw_b,Span_b,SM_Method_b);
error_Smooth_b = (abs(psd_raw_b - psd_Smooth_b)).^2;

psd_Smooth = [psd_Smooth_a; psd_Smooth_b];
error_Smooth = [error_Smooth_a; error_Smooth_b];
