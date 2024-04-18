function filtdata = eqbandpassfir(data,filt_low,filt_high,dt,filt_order)
%	Function to perform firfilter on "data"
%	vector,	and return filtered data vector	"filtdata"
%
%	filter corners:
%         filt_low = low pass
%         filt_high = high pass
%	dt = sampling rate (in sec/sample)
%	filt_order = filter order 
%
samp = 1/dt;		% sampling frequency
%ch1n = ch1*2/samp;
%ch2n = ch2*2/samp;

%disp(['DEBUG: ch1n=',num2str(ch1n),' ch2n=',num2str(ch2n),' samp=',num2str(samp),' dt=',num2str(dt)]);
disp(['DEBUG: filt_low=',num2str(filt_low),' filt_high=',num2str(filt_high),' samp=',num2str(samp),' dt=',num2str(dt)]);
D1 = designfilt('bandpassfir','FilterOrder',filt_order,...
     'CutoffFrequency1',filt_low,'CutoffFrequency2',filt_high,...
     'SampleRate',dt);
filtdata = filter(D1,data);
