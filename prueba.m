
names = dbnames(sl);

for i = 1:length(names)

    ind = regexp(names(i), '^L_', 'match');
        if ~isempty(ind{1})
            plot(sl.(names{i}))
        end
end

