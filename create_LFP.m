function [LFP, powhil_rc] = create_LFP( spikes, ns, pband, sim_length, pre_stim_length )

    SR=1000;% sampling rate
    time=0:1/SR:(sim_length)/1000;
    MUA = zeros(ns, sim_length+1);
    %% add spike events to MUA
    for m=1:size(spikes,1)
        MUA(spikes(m,1),spikes(m,2))=MUA(spikes(m,1),spikes(m,2))+1;
    end
    
    MUA_PSP = zeros(size(MUA,1),size(MUA,2));
    PSP = compute_psp(1.5);
     %% add EPSP to MUA data
    for m1=1:size(MUA,1)
        for m2 = 1:size(MUA,2)
           if(MUA(m1,m2)==1) 
               for l = 1:length(PSP)
                  if(m2+l <= size(MUA,2))
                     MUA_PSP(m1, m2+l) = MUA_PSP(m1, m2+l) + PSP(l);
                  end
               end
           end
        end
    end
    
    %% sum MUA & Hanning filter on LFP
    LFP=sum(MUA_PSP,1);
    filtord=30;
    b=hann(filtord)./filtord;
    LFP=filtfilt(b,1,LFP); 
    
    %% filter LFP by frequency
    if(pband(2) - pband(1) <= 6) % for small bands
        Fstop1 = 0.5; 
        Fstop2 = 101;           % Second Stopband Frequency
        Astop1 = 60;          % First Stopband Attenuation (dB)
        Apass  = 1;           % Passband Ripple (dB)
        Astop2 = 80;
    else % for larger bands
        Fstop1 = 0.00001;
        Fstop2 = 80;           % Second Stopband Frequency
        Astop1 = 40;          % First Stopband Attenuation (dB)
        Apass  = 1;           % Passband Ripple (dB)
        Astop2 = 80;
    end
    Fs = 1000;  % Sampling Frequency
    Fpass1 = pband(1);           % First Passband Frequency
    Fpass2 = pband(2);           % Second Passband Frequency
    match  = 'stopband';  % Band to match exactly
    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                          Astop2, Fs);
    Hd = design(h, 'butter', 'MatchExactly', match);
    LFP=filtfilt(Hd.sosMatrix,Hd.ScaleValues,LFP);
        
    
    %% pre-post stim change
    powhil=abs(hilbert(LFP));
    bl=[0.25 (pre_stim_length/1000)-0.25];
    bl1=find(bl(1)==time);
    bl2=find(bl(2)==time);
    powhil_rc=(powhil-mean(powhil(bl1:bl2),2))./mean(powhil(bl1:bl2),2);
end

