function treatments = loadEntrainTreatments(filename, varargin)

treatdata = h5read(filename,'/Output/Command');
names = h5readatt(filename,'/Output/Command','Channels');
sigtypes = h5readatt(filename, '/Output/Command','SignalTypes');

for i = 1:length(names)
    fldname = genvarname(lower(names{i}));
    treatments.(fldname) = treatdata(:,i);
end

treatments.typeind = treatments.type + 1;
treatments.type = sigtypes(treatments.typeind);


