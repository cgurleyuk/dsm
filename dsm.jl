module dsm
    import DSP

    function snr(y, f)
        n = length(y);
        
        xs = sin.(2*pi*f*(1:n));
        xc = cos.(2*pi*f*(1:n));
        a  = sum(2*xs.*y)/n;
        b  = sum(2*xc.*y)/n;

        i_sig = (a*xs+b*xc);
        i_noi = y - i_sig;
        
        w     = DSP.Windows.kaiser(n,20);
        wf    = sum(w.^2)/n;
        I_sig = fft(i_sig.*w)/n;
        I_noi = fft(i_noi.*w)/n;
        
        P_sig = sum(abs.(I_sig[1:end]).^2);
        P_noi = sum(abs.(I_noi[1:end]).^2);
        
        10*log10(P_sig/P_noi)
    end
end
