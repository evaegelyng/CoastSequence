The folder contains six text files with environmental data related to the sediment and water samples.
The sample IDs are composed as follows: [Sampling_campaign]_[Sampling_cluster]_[Habitat]_[Substrate]_[Sample_replicate],
where the "Sampling_campaign" part is omitted for samples from the first campaign (corresponding to the spring season) and a "2" is included for
samples from the second campaign (corresponding to the autumn season). The "Sampling_cluster" refers to the sampling station and is indicated by a "C"
followed by an integer. For the "Habitat", and "E" indicates eelgrass, while an "R" indicates rocks and an "S" indicates sand bottom. 
The "Substrate" part indicates whether the sample was taken from the sediment ("B" for bottom) or water ("W").
The sample replicate is only included for sediment samples, where the environmental measurements directly correspond to individual sediment samples (replicates 1-3). 
In contrast, measurements of water isotopes, chlorophyll and nutrients were collected from dedicated samples taken in parallel with eDNA samples.

The files are described individually below:

1. CNdelta_iso_sed_both_new.txt
    - Description: Isotope ratios for carbon and nitrogen, respectively, that were measured for sediment samples.
    - Columns: "Sample_ID" (see explanation above), "d14N_15N" (N isotope ratio), "Norm_d14N_15N" (normalized N isotope ratio), "d12C_13C" (C isotope ratio),
    "Norm_d12C_13C" (normalized C isotope ratio), "Remarks" (lab notes; "QL" is quantification limit and "DL" is detection limit), "QC" (quality control; values of "OK" indicate 
    that N and C were both above the quantification limit, while "N_und_QL" indicates that the N content was below the quantification limit etc.) 
2. CNT_sed_both_new.txt
    - Description: Total carbon and total nitrogen concentrations measured for sediment samples, as well as C:N ratio.
    - Columns: "Sample_ID" (see explanation above), "N_per_kg_dry_sample" (N concentration), "C_per_kg_dry_sample" (C concentration), "CN_ratio", (C to N ratio),
    "Remarks" (lab notes; "QL" is quantification limit and "DL" is detection limit), "QC" (quality control; values of "OK" indicate 
    that N and C were both above the quantification limit, while "N_und_QL" indicates that the N content was below the quantification limit etc.) , and whether C and/or N levels were below the quantification limit.
3. dOdH_wat_both:
    - Description: Isotope ratios for hydrogen and oxygen, respectively, that were measured for water samples.
    - Columns: "Sample_ID" (see explanation above), "d18O" (average O isotope ratio across technical replicates), "SD_d18O" (standard deviation for O isotope ratio across technical replicates),
    "d12H" (average H isotope ratio), "SD_d12H" (standard deviation for H isotope ratio)
4. NUT_ctd_wat_both.txt:
    - Description: Nutrient concentrations (all in umol/L) and chlorophyll concentrations measured for water samples, as well as salinity and temperature measurements collected via CTD casts. 
    - Columns: "Sample_ID" (see explanation above), "NO2NO3" (summed nitrite and nitrate concentration), "NO2" (nitrite conc.), "NH3" (ammonia conc.), "PO4" (phosphate conc.), "Si" (silicate conc.),
    "NO3" (nitrate conc.), "Latitude" (GPS-based latitude), "Longitude" (GPS-based longitude), "Chlorophyll" (chlorophyll conc. in yg/L), "Salinity" (salinity ratio using PSS), "Temperature" (degrees celcius), 
    "CTD_depth" (depth in cm of CTD measurements), "Time" (date of sampling).
5. TP_org_inorg_dens_sed_autumn.txt
    - Description: Physicochemical measurements for sediment samples collected in autumn.
    - Columns: "Sample_ID" (see explanation above), "Grain size" (grain size category; 1 = mud, 2= sand, 3 = gravel), "Density_wet" (wet density in g/mL), "Dry_matter_per_WW" (g dry matter per g wet weight),
    "Watercontent" (water content in g/mL), "Organic_content" (g organic matter per g sample), "Inorganic_content" (g inorganic matter per g sample), "TP" (total phosporus in mg/g dry sample).
6. TP_org_inorg_dens_sed_spring.txt
    - Description: Physicochemical measurements for sediment samples collected in spring.
    - Columns: "Sample_ID" (see explanation above), "Grain size" (grain size category; 1 = mud, 2= sand, 3 = gravel), "Density_wet" (wet density in g/mL), "Dry_matter_per_WW" (g dry matter per g wet weight),
    "Watercontent" (water content in g/mL), "Organic_content" (g organic matter per g sample), "Inorganic_content" (g inorganic matter per g sample), "TP" (total phosporus in mg/g dry sample).