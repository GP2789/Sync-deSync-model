function spikes = spike_train( spikesPerS, duration_s, numTrains)

t = [0:1/1000:duration_s/1000-1/1000];
total_spikes = zeros(numTrains, length(t));

for train = 1:numTrains
    vt = rand(size(t));
    total_spikes(train, :) = (spikesPerS*(1/1000)/numTrains) > vt;
end
spikes = sum(total_spikes,1);

end