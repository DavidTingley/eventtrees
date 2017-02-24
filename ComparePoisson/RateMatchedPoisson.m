load('spikes.mat')

numspikes = cellfun(@length,spikes);
numcells = length(numspikes);
% [~,sortrate] = sort(numspikes);

dt = 5; %s
overlap = 1;
window = dt.*overlap;
[spikemat,t,spindices] = SpktToSpkmat(spikes, [], dt,overlap);
ratemat = spikemat./window;

%%
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

    %%
    figure
    imagesc(t,1:numcells,log10(ratemat)')
    axis xy
    hold on
    plot(poisspikes_t,poisspikes_cell,'r.')
    ColorbarWithAxis([-1,1],'Rate (logHz)')
    xlim([0 20])
    xlabel('t (s)');ylabel('Cell (sorted by mean rate)')
    