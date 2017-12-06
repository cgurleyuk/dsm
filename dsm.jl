module dsm
    export Mod1, Mod2
    export snr, sim_mod, sim_snr

    import DSP

    abstract type Mod end

    type Mod1 <: Mod
        a_1::Float64
        b_1::Float64
    end

    type Mod2 <: Mod
        a_1::Float64
        a_2::Float64
        b_1::Float64
        b_2::Float64
        c_1::Float64
        c_2::Float64
        d_13::Float64
    end

    function uniform_quantizer(n::Int64)
        return collect(linspace(-1, 1, n+1)[2:end-1]);
    end

    # TODO: implement random quantizer
    function random_quantizer(n::Int64)
        return collect(linspace(-1, 1, n+1)[2:end-1]);
    end

    function quantize(in::Float64, q::Array{Float64, 1})
        return 2*(searchsorted(q, in).stop/length(q)-.5);
    end

    function quantize(in::Array{Float64, 1}, q::Array{Float64, 1})
        nl = length(q);
        r = zeros(size(in));
        for i = 1:length(in)
            r[i] = 2*((searchsorted(q, in[i]).stop)/nl-.5);
        end

        return r
    end

    function sim_mod(md::Mod2, in::Array{Float64, 1})
        szin = size(in)
        int1 = zeros(szin);
        int2 = zeros(szin);
        q    = zeros(szin);
        q[1] = 1;

        for i = 2:length(in)
            int1[i] = int1[i-1] + md.b_1*in[i-1] - md.a_1*q[i-1];
            int2[i] = int2[i-1] + md.c_1*int1[i] + md.b_2*in[i-1] - md.a_2*q[i-1];
            q[i] = ((md.c_2*int2[i] + md.d_13*int1[i]) > 0)*2-1;
        end

        return (q, int1, int2)
    end

    function sim_mod(md::Mod1, in::Array{Float64, 1})
        szin = size(in);
        int1 = zeros(szin);
        y    = zeros(szin);
        q    = uniform_quantizer(2);

        y[1] = 1;

        for i = 2:length(in)
            int1[i] = int1[i-1] + md.b_1*in[i-1] - md.a_1*y[i-1];
            y[i] = quantize(int1[i], q);
        end

        return (y, int1)
    end

    function sim_snr(md::Mod, a)
        snr = zeros(size(a));

        n_fft = 2^16;
        n_i   = 2^5;
        f_i   = 1/n_fft*57;
        f_bw  = 1/(2*256);
        t     = collect(0:(n_fft+n_i-1));

        for i = 1:length(a)
            inp    = a[i]*sin.(2*pi*f_i*t);
            out    = dsm.sim_mod(md, inp);
            
            q      = out[1][end-n_fft+1:end];
            snr[i] = dsm.snr(q, f_i, f_bw, 1);
        end

        return snr
    end

#    function snr_iq(y, f_i, f_bw)
#        n   = length(y);
#        nbw = trunc(Int64, n*f_bw);
#
#        xs = sin.(2*pi*f_i*(1:n));
#        xc = cos.(2*pi*f_i*(1:n));
#        a  = sum(2*xs.*y)/n;
#        b  = sum(2*xc.*y)/n;
#
#        i_sig = (a*xs+b*xc);
#        i_noi = y - i_sig;
#
#        w     = DSP.Windows.hanning(n);
#        wf    = sum(w.^2)/n;
#        wfa   = sum(w)/n;
#
#        I_sig = fft(i_sig.*w)/n
#        I_noi = fft(i_noi.*w)/n;
#
#        P_sig = sum(abs.(I_sig[1:nbw]).^2);
#        P_noi = sum(abs.(I_noi[1:nbw]).^2);
#
#        10*log10(P_sig/P_noi)
#    end

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

