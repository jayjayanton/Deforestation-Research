# Input datasets

I haven't committed the original WWF datasets. The portable Colombia script expects a CSV with `Deforestation` coded 0/1, numeric predictors, and an `Elevation` column. The Brazil script expects a CSV with `DEFORESTATION` coded 0/1, numeric predictors, and `ELEVATION_SRTM90_AMAZONIA`.

The original Brazil workbook was `Amazonia_2004_2017_03_24_sel_coord.xlsx`. Convert an authorized copy to CSV before running the portable script. For the Colombia workflow, the original input was `Colombia_Analysis_03_07_2021_No_IDs.csv`.

Review predictor units and label definitions before running. Use `--predictors` to explicitly choose predictors, and `--distances-in-metres true` only when appropriate. The source data are not publicly distributed here.
