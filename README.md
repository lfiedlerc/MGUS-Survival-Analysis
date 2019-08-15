# Estimating Survivorship for Monoclonal Gammopathy of Undetermined Significance

## Objective
The objective of this project is to model the progression of individuals with monoclonal gammopathy of undetermined significance to a plasma cell malignancy, and to identify the significant predictors of a plasma cell malignancy.

MGUS is a condition where there is an abnormal protein in the blood (Monoclonal gammopathy, 2017). The monoclonal immunoglobulin (Ig), or M protein, is produced by the plasma cells in the bone marrow. In MGUS, the M protein can accumulate to such levels that it inhibits healthy cells and can lead to tissue damage. Although the condition is generally asymptomatic and very seldom problematic, MGUS can progress to more serious disorders such as blood cancer. In fact, the study done by van de Donk et al. (2014) found that MGUS will commonly precede multiple myeloma, a cancer of the plasma cells.

## Methods
Patient records collected between 1960 and 1994 for 1384 individuals with monoclonal gammopathy of undetermined significance were analyzed using R-studio. Kaplan-Meier was used to obtain a crude survivorship; subsequently a Cox proportional hazard model was used to identify the significant predictors of plasma cell malignancy.

### Dataset
The dataset that will be used for this project has been donated courtesy of Dr. Robert Kyle of the Mayo Clinic (Kyle et al., 2002). It contains 1341 records of de-identified data describing the natural history of patients with MGUS.

Each observation contains the following 10 variables:

1. ID: A numeric identifier for the patient
2. Age: The age of diagnosis in years
3. Sex: The gender of the patient
4. HGB: Hemoglobin values
5. Creat: Serum creatinine values
6. MSpike: Size of the monoclonal serum spike
7. Ptime: Time until progression to a plasma cell malignancy (PCM) or last contact, in months
8. Pstat: Binary indicator of a PCM event
9. Futime: Time until death or last contact, in months
10. Death: Binary indicator of a death event

## Results
Median time before development of a plasma cell malignancy was 31 years. Hgb and monoclonal serum spike were found to be significant predictors of progression to plasma cell malignancy. A hazard ratio of 2.48 per unit increase was estimated for serum monoclonal concentration levels, as well as an 11% reduction in hazard per unit increase of Hemoglobin levels.

## Conclusions
High survival probabilities even at longer time points reflect the low prevalence of progression to a plasma cell malignancy. This study is limited by the use of potentially outdated data and by the low number of events, resulting in poor statistical power to detect possible differences in survivorship between genders. Future work should consider the most recent criteria for diagnosing monoclonal gammopathy of undetermined significance as well as the inclusion of race and family history as a covariate.

## References
-Monoclonal gammopathy of undetermined significance (MGUS). (2017, July 29). Retrieved March 03,
	2018, from https://www.mayoclinic.org/diseases-conditions/mgus/symptoms-causes/syc-20352362
-Kyle, R. A., Therneau, T. M., Rajkumar, S. V., Offord, J. R., Larson, D. R., Plevak, M. F., & Melton III,
	L. J. (2002). A long-term study of prognosis in monoclonal gammopathy of undetermined significance. New England Journal of Medicine, 346(8), 564-569.
-van de Donk, N. W., Palumbo, A., Johnsen, H. E., Engelhardt, M., Gay, F., Gregersen, H., ... & Musto, P.
	(2014). The clinical relevance and management of monoclonal gammopathy of undetermined significance and related disorders: recommendations from the European Myeloma Network. Haematologica, 99(6), 984-996.
