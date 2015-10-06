
b__vec_ = [.002 ]; % 5 points
rImit__ = 0.1;

for icd = 1:length(b__vec_)
    bHarvest__ = b__vec_(icd);

    full_model_v3_1_LI_Het_50;

    plot_stuff;
end

