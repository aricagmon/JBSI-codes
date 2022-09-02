function syncData = jbsi_synchrony(a, b, sw)
% jbsi_synchrony  JBSI synchrony analysis of Aric Agmon
% 
%    syncData = jbsi_synchrony(a, b, sw) computes the Jitter-Based 
%    Synchrony Index between spike trains a and b for the synchrony 
%    window sw. 
% 
%    a and b are vectors of spike times, in seconds. 
% 
%    sw is the synchrony window width, in ms. sw may be a scalar or vector. 
%    If sw is not input, then the following values are used: 
%       sw = sqrt(2).^(-3:14), or 0.35 to 128 ms.
%
%    syncData is a struct holding the results.
% 
%    The algorithm used to calculate the JBSI was created by Aric Agmon of
%    West Virginia University and may be found in 
% 
%          Agmon (2012). A novel, jitter-based method for detecting
%          and measuring spike synchrony and quantifying temporal firing 
%          precision. Neural Systems and Circuits, 2:5.
% 
%    The algorithm was originally coded in MathCad by Aric Agmon. 
%    The algorithm was translated to Matlab by
% 
%    Craig Atencio
%    9/16/2012


error(nargchk(2,3,nargin));

if ( nargin == 2 ) % define synchrony window since it wasn't provided by user
   n = -3:14; % use for sqrt(2).^(n)
   sw = sqrt(2) .^(n); % allows for a variety of synchrony windows

%    n = -2:7; % use for 2.^(n)
%    sw = 2 .^(n); % allows for a variety of synchrony windows

   jw = 2*sw;
else % the user supplied the synchrony window
   jw = 2*sw;
end



%-------------------------------------------------------------------
%-------------------------------------------------------------------
% I've left this in to match the MathCad code. This can be modifed
% for a general version of the code. 
%-------------------------------------------------------------------

% start and finish times. Don't use end since that's a reserved word
start0 = 0;
finish0 = Inf; % can make this 100, or 260 to include less data


% Make sure spikes are in ascending order
a = sort(a);
b = sort(b);


% Only keep spike times > 0
a = a(a>0);
b = b(b>0);


% Get spike times that are between the start and finish times
a = a( a>start0 & a<finish0 );
b = b( b>start0 & b<finish0 );

%-------------------------------------------------------------------
%-------------------------------------------------------------------



% The Reference spike train has fewer spikes.
% The Target has more spikes.
%
% Assign to variables Reference and Target so that the names are more
% descriptive.

if ( length(a) <= length(b) )
   Reference = a;
   Target = b;
else
   Reference = b;
   Target = a;
end

fprintf('\n');
fprintf('Spikes in Reference Train = %.0f\n', length(Reference));
fprintf('Spikes in Target Train = %.0f\n', length(Target));


% Check to make sure no spikes are separated by less than 2 ms
% --------------------------------------------------------------
% dRef = 1000*diff(Reference);
% indRef = find(dRef < 2);
% checkReference = length(indRef);
% 
% if ( checkReference )
%    error('Error in Reference spike train.');
% end
% 
% 
% dTar = 1000*diff(Target);
% indTar = find(dTar < 2);
% checkTarget = length(indTar);
% 
% if ( checkTarget )
%    error('Error in Target spike train.');
% end



sync_frac = zeros(1,length(sw));
aveNc = zeros(1,length(sw));
varNc = zeros(1,length(sw));
sdNc = zeros(1,length(sw));
zscore = zeros(1,length(sw));
beta = zeros(1,length(sw));
jbsi = zeros(1,length(sw));

fprintf('\n');
fprintf('SW\t\tJW\n');

for i = 1:length(sw)

   fprintf('%.3f\t%.3f ms\n', sw(i), jw(i) );

   % Shift Target spike train to calculate synchrony for non-zero lags
   % --------------------------------------------------------------
   shift = 0;
   Target = shift_target_spikes(Target, shift);


   % Calculate the non-overlapping synchrony windows
   % --------------------------------------------------------------
   seg = synchrony_windows(Target, sw(i));


   % For every spike in the Reference train, calculate the probability 
   % that after a jitter JW, the spike will be synchronous with a segment
   % in the Target train. This gives the probability for random
   % synchronous spiking
   % --------------------------------------------------------------
   psi = synchrony_probability(Reference, seg, jw(i));


   % Calculate the number of spikes in the synchrony window
   % --------------------------------------------------------------
   Nc(i) = synchronous_spikes(Reference, Target, sw(i)); % coincidence count, Nc


   % The following computes the JBSI without using the explicit PDF:
   % --------------------------------------------------------------
   sync_frac(i) = Nc(i) / length(Reference); % Index


   % Use calculated PDF values to obtain the JBSI
   %---------------------------------------------------------------
   aveNc(i) =  sum(psi); % average count, <Nc>, expected by chance
   varNc(i) = sum( psi .* (1-psi) );
   sdNc(i) = sqrt(varNc(i));
   zscore(i) = (Nc(i) - aveNc(i)) / sdNc(i); % z-score


   % Get beta value
   if ( jw(i) / sw(i) <= 2 )
      beta(i) = 2;
   else
      beta(i) = jw(i) / (jw(i) - sw(i));
   end


   % Calculate jbsi using pdf: jbsi = beta * (Nc - <Nc>) / n1
   %---------------------------------------------------------------
   jbsi(i) = beta(i) * (Nc(i) - aveNc(i)) / length(Reference); % jitter based sensitivity index


%    % Uncomment if you want to print out the values while processing
%    % --------------------------------------------------------------
%    fprintf('index = %.3f\n', sync_frac(i));
%    fprintf('Nc = %.0f\n', Nc(i));
%    fprintf('<Nc> = %.3f\n', aveNc(i));
%    fprintf('Var(Nc) = %.3f\n', varNc(i));
%    fprintf('SD(Nc) = %.3f\n', sdNc(i));
%    fprintf('z score = %.3f\n', zscore(i));
%    fprintf('beta = %.3f\n', beta(i));
%    fprintf('jbsi = %.3f\n', jbsi(i));

end % (for i)


maxISI = 100; % maximum delay between spikes
cISI = get_cISI(Reference, Target, maxISI);

binwidth = 2; % ms, for the correlogram
[delay, ccraw, ccflat, cc] = get_crosscorr(cISI, maxISI, binwidth, Reference, Target);


% Assign variables to output struct
syncData.sw = sw;
syncData.jw = jw;
syncData.n1 = length(Reference);
syncData.n2 = length(Target);
syncData.minspike = min([min(Reference) min(Target)]);
syncData.maxspike = max([max(Reference) max(Target)]);
syncData.sync_frac = sync_frac;
syncData.Nc = Nc;
syncData.aveNc = aveNc;
syncData.varNc = varNc;
syncData.sdNc = sdNc;
syncData.zscore = zscore;
syncData.beta = beta;
syncData.jbsi = jbsi;
syncData.maxISI = maxISI;
syncData.cISI = cISI;
syncData.binwidth = binwidth;
syncData.delay = delay;
syncData.ccraw = ccraw;
syncData.ccflat = ccflat;
syncData.cc = cc;


% Plot if more than one jbsi value calculated
if ( length(sw) > 1 )
   plot_cc_jbsi(syncData);
end


return;



%-------------------------------------------------------------------
% Function Definitions
%-------------------------------------------------------------------

function syn = synchronous_spikes(Reference, Target, SW)
% synchronous_spikes  Count of observed synchronous spikes
% 
%    syn = synchronous_spikes(Reference, Target, SW) 
% 
%    Counts the number of observed synchronies for the reference spike
%    train Reference and the target spike train Target, where 
%    length(Reference) < length(Target. Each vector of spike times is 
%    in seconds. Spikes are determined to be synchronous if they fall 
%    within the synchrony window, SW, in ms.
% 
%    This algorithm counts only the first Target spike that is synchronous 
%    with any Reference spike, so it produces correct results even if 
%    SW is larger than 0.5*min(ISI_Target) and synchrony windows overlap. 
%    Therefore, each spike in the reference train can add no more 
%    than 1 to the synchrony count.
% 
%    Reference : vector of spike times, in seconds
% 
%    Target : vector of spike times, in seconds
% 
%    SW : synchrony window, in ms
%
%    syn : count of synchronous spikes. An integer.


syn = 0;

for i = 1:length(Reference)

   found = 0;
   j = 1;

   while ( j <= length(Target) && ~found )

      if ( 1000 * abs( Reference(i) - Target(j) ) <= SW )
         syn = syn + 1;
         found = 1;
      end

      j = j + 1;

   end % (while j)

end % (for i)

return;



function psi = synchrony_probability(Reference, seg, JW)
% synchrony_probability  Probability of synchrony for a jittered Refernce spike 
% 
%    psi = synchrony_probability(Reference, seg, JW)
% 
%    Calculates the probability psi(i) that after a jitter JW, spike i 
%    in the reference train will be synchronous with a segment in the 
%    Target train.
% 
%    Reference : vector of spike times, in seconds, for Reference spike train
% 
%    seg : matrix of non-overlapping synchrony windows
% 
%    JW : jitter window, in ms

error(nargchk(3,3,nargin));


psi = zeros( 1, length(Reference) ); % probability value for each Reference spike

for i = 1:length(Reference) % go through reference train spikes

   for h = 1:size(seg,1) % go through all synchrony window segments

      if ( (1000 * seg(h,1) - JW) <= 1000*Reference(i) && ...
            1000*Reference(i) <= min([ 1000*seg(h,1)+JW   1000*seg(h,2)-JW ]) )

         psi(i) = psi(i) + (JW - 1000*(seg(h,1)-Reference(i) ) ) / (2*JW);


      elseif ( 1000*Reference(i) > min([ 1000*seg(h,1)+JW  1000*seg(h,2)-JW ]) && ... 
               1000*Reference(i) < max([ 1000*seg(h,1)+JW  1000*seg(h,2)-JW ]) )

         psi(i) = psi(i) + min([( 1000*seg(h,2) - 1000*seg(h,1) )/(2*JW)  1]);

      elseif ( 1000*Reference(i) > max( [ 1000*seg(h,1)+JW  1000*seg(h,2)-JW ] ) && ...
               1000*Reference(i) < (1000*seg(h,2) + JW ) )

         psi(i) = psi(i) + ( JW + 1000 * (seg(h,2) - Reference(i) ) ) / (2*JW);

      else

         psi(i) = psi(i) + 0; % Don't really need this; include for symmetry
                              % and to follow paper and original code
      end

   end % (for h)

end % (for i)

return;



function seg = synchrony_windows(Target, SW)
% synchrony_windows  Set of non-overlapping synchrony windows of Target spike train
% 
%    seg = synchrony_windows(Target, SW)
% 
%    Takes the spikes in the target spike train Target and calculates
%    the set of all synchrony windows having the width SW. 
% 
%    The windows do not overlap, and are later used to calculate the 
%    expected number of chance coincidences between the Target train and 
%    the Reference Train.
% 
%    Target : Target spike train, a vector, with spikes in units of seconds.
% 
%    SW : synchrony window, in ms

error(nargchk(2,2,nargin));


% Calculate synchrony windows
% --------------------------------------------------------------------
% Note: If SW is larger than 1/2 of min(ISI_Target), synchrony windows will 
% overlap, and u will be smaller than r.

g = 1; % index into row of union of synchrony windows
j = 1; % index into spike train Target
done = 0;

while ( j <= length(Target)-1 && ~done)

   seg(g,1) = (1000 * Target(j) - SW) / 1000;

   if ( j == length(Target)-1 )
      done = 1;
   end

   while ( 1000*(Target(j+1) - Target(j)) < 2*SW && ~done)

      j = j + 1;

      if ( j == length(Target)-1 )
         done = 1;
      end

   end % (while)

   seg(g,2) = ( 1000*Target(j) + SW ) / 1000;

   g = g + 1;
   j = j + 1;

end % (while)


return;



function Target = shift_target_spikes(Target, shift)
% shift_target_spikes Shift spikes times in spike train
% 
%    Target = shift_target_spikes(Target, shift)
% 
%    Target : vector of spikes times, in seconds.
%    shift : amount of time, in ms, that Target is to be shifted. 

error(nargchk(2,2,nargin));


% Shift Target spike train by shift ms

Target_shifted = (1000*Target + shift) / 1000; % shift Target

Target = Target_shifted; % reassign Target

return;



function cISI = get_cISI(Reference, Target, maxISI)
% get_cISI  Cross interspike interval function
% 
%    cISI = get_cISI(Reference, Target, maxISI)
%
%    Reference : reference spike times, in seconds.
%    Target : target spike times, in seconds.
%    maxISI : maximum ISI to process
% 
%    cISI : vector of interspike times

error(nargchk(3,3,nargin));


k = 1;

for i = 1:length(Reference)

   for j = 1:length(Target) 

      dif = 1000 * ( Reference(i) - Target(j) );

      if ( dif > -maxISI && dif < maxISI )
         cISI(k) = dif;
         k = k + 1;
      end

   end % ( for j )

end % (for m)

return;



function [delay, ccraw, ccflat, cc] = get_crosscorr(cISI, maxISI, bw, Reference, Target)
% get_crosscorr Correlogram from Interspike Intervals
% 
%    [delay, ccraw, ccflat, cc] = get_crosscorr(cISI, maxISI, binwidth, Reference, Target)
%
%    Takes the crossISI data in cISI and computes the correlogram for maximum
%    delay maxISI and binning width bw, in ms. Reference and Target are
%    the spike trains, in units of seconds.
%
%    cISI is obtained from the function get_cISI.

error(nargchk(5,5,nargin));


% Delays for correlograms
binpos = 0:bw:maxISI;
delay = [-binpos(end:-1:2) binpos];

% Correlogram without normalization based on inter-spike intervals
ccraw = hist(cISI, delay);


% Total duration of spike train recording in milliseconds
T = 1000 * ( max([max(Reference) max(Target)]) - min([min(Reference) min(Target)]) );

% Normalizing factor from correlogram
nf = (T - abs( delay ) ) / T;

% Correlogram normalized
ccflat = ccraw ./ nf;

% Correlation coefficient
cc = ccflat / length(Reference);

return;



function plot_cc_jbsi(syncData)
% plot_cc_jbsi  Graph Correlation coefficient and JBSI
% 
%    plot_cc_jbsi(syncData) takes data in the struct syncData and plots
%    the correlation coefficient versus delay. It then plots the JBSI
%    versus synchrony window, or precision.

error(nargchk(1,1,nargin));


% Get correlation parameters
delay = syncData.delay;
cc = syncData.cc;
maxISI = syncData.maxISI;
binwidth = syncData.binwidth;

% Get jbsi parameters
sw = syncData.sw;
zscore = syncData.zscore;
jbsi = syncData.jbsi;


% Plot the results
close all;

figure;

subplot(2,1,1);
hb = bar(delay, cc);
set(hb,'facecolor', 0*ones(1,3));
set(hb,'edgecolor', 0*ones(1,3));
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlim([-maxISI-binwidth maxISI+binwidth]);
ylim([0 1.1*max(cc)]);
box off;
xlabel('Delay (ms)');
ylabel('CC');
title(sprintf('Binwidth = %.2f ms', binwidth));


subplot(2,1,2);
plot(sw, jbsi, 'ks-', 'markerfacecolor', 'k', 'markersize', 3);
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlim([min(sw)-0.05*max(sw) max(sw)+0.05*max(sw)]);
box off;
xlabel('Precision (ms)');
ylabel('JBSI');

set(gcf,'position', [680   266   660   712]);

return;
