function [iwant,nval] = sub_percentile(x,Pwant);
%
%  For our Nfound data we can not use STD statistical values because the data set is not guassian, they have long tails!
%  Instead we will use percentile looking for the percentile above Pwant 
%
Nx = length(x);
%
%  If data contains negative values shift .... 
%
if min(x)<0
   disp('Caution: Data contains negative values!');
   keyboard
end
x_round = round(x);
dx = max(x_round)-min(x_round)
dstep = dx/100;
%
%  Compute the percentiles
%
maxV = max(x_round);
disp('Double Check the below!');
for i=1:1:maxV+1
%for i=1:1:100
    value(i) = i-1;
    %value(i) = maxV-(i-1)*dstep;
    nn(i) = length(find(x_round==value(i)));
    %value(i) = maxV-(i-1)*dstep;
    %nn(i) = length(find(x>value(i)));
end

perC = nn./Nx;

cumsumperC = cumsum(perC);
nval = min(find(cumsumperC>=Pwant/100));
%keyboard
disp(['DEBUG: nval=',num2str(nval),' max(x)=',num2str(max(x))]);
%keyboard
iwant = find(x>=nval);

iplot = 1;
if iplot==1
   figure(33)
   clf
   subplot(2,1,1)
   plot(x,'g.')
   hold on
   plot(iwant,x(iwant),'mo')
   legend('All',[num2str(Pwant),'%'])
   subplot(2,1,2)
   plot(cumsumperC,'r.');
   hold on
   yy = get(gca,'Ylim');
   plot([nval nval],yy,'b-');
   title(['nval=',num2str(nval)]);
   %disp('Hit any key to continue... ');
   %pause
end
%pjunk = ['plot_',num2str(Pwant),'.png'];
%print(pjunk,'-dpng')
