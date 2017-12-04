module dsm
    export snr, snrx

    import DSP

    function snr_iq(y, f_i, f_bw)
        n   = length(y);
        nbw = trunc(Int64, n*f_bw);

        xs = sin.(2*pi*f_i*(1:n));
        xc = cos.(2*pi*f_i*(1:n));
        a  = sum(2*xs.*y)/n;
        b  = sum(2*xc.*y)/n;

        i_sig = (a*xs+b*xc);
        i_noi = y - i_sig;

        w     = DSP.Windows.hanning(n);
        wf    = sum(w.^2)/n;
        wfa   = sum(w)/n;

        I_sig = fft(i_sig.*w)/n;
        I_noi = fft(i_noi.*w)/n;

        P_sig = sum(abs.(I_sig[1:nbw]).^2);
        P_noi = sum(abs.(I_noi[1:nbw]).^2);

        10*log10(P_sig/P_noi)
    end

    function snr_kaiser(y, f_i, f_bw, f_s)
        n = length(y);

        sig_bin  = f_i/f_s*n;
        bw_bin = f_bw/f_s*n;
        sig_bins = (sig_bin-14:1:sig_bin+16);
        noi_bins = vcat(1:sig_bin-15, sig_bin+17:bw_bin);

        w = DSP.Windows.kaiser(n, 20);
        wf    = sum(w.^2)/n;
        Y = fft(y.*w)/n;
        P_sig = sum(abs.(Y[trunc.(Int64, sig_bins)]).^2);
        P_noi = sum(abs.(Y[trunc.(Int64, noi_bins)]).^2);

        10*log10(P_sig/P_noi)
    end

    function snr(y, f_i, f_bw, f_s)
        # snr_iq(y, f_i/f_s, f_bw/f_s);
        snr_kaiser(y, f_i, f_bw, f_s);
    end

end

