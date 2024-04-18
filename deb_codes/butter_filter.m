function [data_filt] = butter_filter(data_in, lo_corner, hi_corner, dt, norder, ftype ) 
%
% Purpose: Filter a data trace using a butterworth filter.
% The parameters lo_corner & hi_corner denote the frequency band
% of interest. The time step is 'dt' and 'norder' is the  
% order of the filter.  'ftype' specifies the type of filter, 
% possibilities include:
%
%     high: high pass 
%     low: low pass
%     bandp: bandpass
% 
% All input frequencies are in Hz.  Inputs to Matlab routines must
%    be scaled by the Nyquist frequency 
%
% Norder is half the actual filter order.  
%
% Fails if lo_corner > hi_corner or if ftype is not set.
%
% Example: butter_filter(res(6).wholewf,0.2, 50,res(6).dt,6,'high');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear A B
nyf = 1/(2*dt);
%
%  The butter function requires hi_corner not be equal to or
%  larger than the Nyquist frequency
%
if (hi_corner >= nyf)
   hi_corner = 0.99999*nyf;
end


if strncmp(lower(ftype), 'bandp', length(ftype))

     if (lo_corner > hi_corner)
        warning('Error in butter_filter function: lo_corner > hi_corner')
        pause on
        return
     end
      [B,A]=butter(norder, [lo_corner hi_corner]/nyf );

elseif strncmp(lower(ftype), 'high', length(ftype))

      [B,A]=butter(norder,lo_corner./nyf,'high');

elseif strncmp(lower(ftype), 'low', length(ftype))
 
      [B,A]=butter(norder, hi_corner/nyf);

else

    disp(['Unknown filter type: ' ftype '  No filter applied'])
    data_filt = data_in ;
    return

end

%
%  Now do the final filter
%
data_filt=filtfilt(B,A,data_in);

return

%
%  Debug plot stuff
%
[junk_spec1,junk,junk_f1] = fftsac(data_in,dt);
[junk_spec2,junk,junk_f2] = fftsac(data_filt,dt);

figure(12)
clf
subplot(2,1,1)
plot(dt.*[0:length(data_in)-1],data_in,'r-')
xlabel('Time (seconds)')
hold on
plot(dt.*[0:length(data_in)-1],data_filt,'c-')
legend('Original','Filtered')
title('Debug Plot')
axis tight
hold off

subplot(2,1,2)
spmax = 1.1 * max(junk_spec2) ;
spmin = spmax/1e3 ;
plot(junk_f1(2:length(junk_f1)),junk_spec1(2:length(junk_spec1)),'r', ...
    junk_f2(2:length(junk_spec2)), junk_spec2(2:length(junk_spec2)),'c')
set(gca,'Ylim',[spmin spmax])
set(gca,'YScale', 'log', 'XScale', 'log')
xlabel('Freq (Hz)')
axis tight
hold off
% keyboard
