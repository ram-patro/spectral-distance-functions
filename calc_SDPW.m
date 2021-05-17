function RSDP=calc_SDPW(a,b,D,type)
RSDP=max([spectral_distance_all(a,D,type)/spectral_distance_all(b,D,type) spectral_distance_all(b,D,type)/spectral_distance_all(a,D,type)]);
end