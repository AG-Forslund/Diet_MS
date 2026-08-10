QCRLSC_CV = 30
Span = 0.75
Rule = 0

write_csv(result_file, here("temp", "temp_feature.csv"))

write_csv(sample_file, here("temp", "temp_sample.csv"))

shiftCor(samPeno = here("temp", "temp_sample.csv"), samFile = here("temp", "temp_feature.csv"), MLmethod = "QCRLSC", QCspan = Span, Frule = Rule, imputeM = "min", coCV = QCRLSC_CV)

result_file_norm <- read.csv(here("Scripts", "statTarget", "shiftCor", "After_shiftCor", "shift_all_cor.csv"))