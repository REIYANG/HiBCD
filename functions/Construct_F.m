function F = Construct_F(HS_spec,MS_spec)

F = zeros(size(MS_spec,1),length(HS_spec));
for count = 1:size(MS_spec,1)
    count_start = find(HS_spec>MS_spec(count,1),1);
    count_end = find(HS_spec<MS_spec(count,2),1,'last');
    F(count,count_start:count_end) = 1/(count_end+1-count_start);
end