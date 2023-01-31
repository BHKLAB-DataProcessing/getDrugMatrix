import pandas as pd
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )

data = pd.read_csv(
    'https://orcestradata.blob.core.windows.net/toxico/DrugMatrix_array_samples.txt')
files = data['x'].values.tolist()

prefix = config["prefix"]

rule get_tset:
    input:
        prefix + 'processed/eset_DM.rds',
        prefix + "data/s_Hepatocyte.csv",
        prefix + "data/curationDrug.rds",
        prefix + "data/doselevel_mapping.rds"
    output:
        prefix + "drugMatrix.rds"
    shell:
        """
        Rscript {prefix}scripts/getDM.R {prefix}
        """

rule process_drug_matrix:
    input:
        expand(prefix + "download/{file}", file=files)
    output:
        prefix + 'processed/eset_DM.rds'
    shell:
        """
        Rscript {prefix}scripts/processDrugMatrixArray.R {prefix}download {prefix}processed
        """

rule download_drug_matrix:
    output:
        expand(prefix + "download/{file}", file=files)
    shell:
        """
        Rscript {prefix}scripts/downloadDrugMatrixArray.R
        """
