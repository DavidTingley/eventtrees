function [poissonspiketimes,spikeplotindices] = RateMatchedPoissonSpikes(spiketimes,dt_rate)
%[poissonspiketimes] = RateMatchedPoissonSpikes(spiketimes,dt_rate)
%simulates poisson spikes with comparable time-varying rate to the spikes
%in spiketimes, with rate time resolution of dt_rate
%
%INPUT
%   spiketimes  {Ncells} cell array of [Nspikes] spiketimes
%   dt_rate     rate time resoluion
%
%OUTPUT
%   poissonspiketimes   poisson simulated spike times, in similar format to
%                       spiketimes
%   spikeplotindices    [Nspikes x 2] vector of plot indicies - first
%                       column is spike times, second is cell number
%
%Aug 2016
%DLevenstein
%% Calculate the Rate Matrix
dt = 1; %s
overlap = 1;
window = dt.*overlap;
[spikemat,t,spindices] = SpktToSpkmat(spikes(sortrate), [], dt,overlap);
ratemat = spikemat./window;

%% Simulate Poisson Spikes for each rate time window
dt_poisson = 0.001;
poisspikes_t = [];
poisspikes_cell = [];

for ww = 1:length(t);
    if mod(ww,1000) == 1
        display(['Window: ',num2str(ww),' of ',num2str(length(t))])
    end
s =  PoissonRateSpikeBins(ratemat(ww,:),dt_poisson,window./dt_poisson);
[s_t,s_cell] = find(s);
s_t = s_t.*dt_poisson+ww-1;

poisspikes_t = [poisspikes_t;s_t];
poisspikes_cell = [poisspikes_cell;s_cell];
end

%% Convert times and cells vectors to cell array format

for cc = 1:numcells
    poissonspiketimes{cc} = poisspikes_t(poisspikes_cell==cc);
end
spikeplotindices = [poisspikes_t,poisspikes_cell];
