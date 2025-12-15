Data used:
- Phenological records from 2019-2024. Only records with data for flowering onset and GSL. The number of data points per species was at least 10.


- Mean annual temperature per garden and year extracted from CRU https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.09/ge/
Flowering Onset: Mean annual temperature + three months before the flowering onset. If flowering was in the second half, I also inlcuded the corresponding month
GSL: Mean annual temperature

Data Cleaning:
I performed a temperature-phenology regression on all data points for a species and eliminated possible outliers using Cook's distance.


Bootstrapping: 
I resampled 1000 times 15 datapoints from the data, without any corrections for specific gardens or something similiar. 
I extracted slopes from the most simplified model lm(Flowering_Onset7GSL ~ Temperature). Again without any implementation of correcetion regarding random effects (i.e. garden or year)

Data summary:
I calculated the mean slope per species and included the lower and upper 95% confindence intervall as well as the variance. I also added per Species the number of Gardens, the min., max. and mean temperature and the number of datapoints.



Files:
Slopes_Bootstraped1000times_all_Flowering
Slopes_Bootstraped1000times_all_GSL

=> Slopes of the 1000 bootstrapping iterations for all Gardens

Overview_Slopes_Flowering
Overview_Slopes_GSL

=> Overview per species including mean slopes
