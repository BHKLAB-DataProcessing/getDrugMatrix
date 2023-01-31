import pandas as pd
data = pd.read_csv(
    'https://orcestradata.blob.core.windows.net/toxico/DrugMatrix_array_samples.txt')
files = data['x'].values.tolist()
