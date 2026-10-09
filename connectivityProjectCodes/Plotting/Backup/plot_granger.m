data = load("/media/venkatesh/3cb61805-2382-46a6-8b40-24cab80d1e2d/ProjectDhyaan/BK1/savedata_2026-10-08_17-43-56/019CKa/EO1_ep_v8_granger.mat");

seed_elec = 1;
data_pre = zeros(64);
for compar_elec = 1:64
data_pre = data_pre + data.connPre(seed_elec,compar_elec,:);
end
data_pre = data_pre/64;
figure;
plot(data_pre);