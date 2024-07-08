import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
import scipy.sparse as sp
import glob
from tqdm import tqdm
import gc  # 导入垃圾回收模块


def load_macaque_data(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0)
    data_matrix = df.pivot_table(
        index='cell_id', columns='gene', values='expr', aggfunc='sum', fill_value=0
    )
    obs_data = (
        df[['chip', 'cell_id', 'gene_area']]
        .drop_duplicates('cell_id')
        .set_index('cell_id')
    )
    spatial_data = (
        df[['cell_id', 'x', 'y', 'rx', 'ry']]
        .drop_duplicates('cell_id')
        .set_index('cell_id')
    )
    obsm_data = {'spatial': spatial_data[['x', 'y', 'rx', 'ry']].values}

    adata = sc.AnnData(
        X=sp.csr_matrix(data_matrix),  # 细胞-基因表达矩阵
        obs=obs_data,  # 观测值矩阵
        var=pd.DataFrame(index=data_matrix.columns),  # 变量矩阵
        obsm=obsm_data,  # 观测值的多维数组
    )

    del df, data_matrix, obs_data, spatial_data, obsm_data
    gc.collect()
    return adata


macaque1 = glob.glob(
    '/home/hanyuji/Data/ST_data/macaque_cortex/ST_result/20230503-macaque1-sct/*macaque1*'
)

h5ad_dic = '/home/hanyuji/Data/ST_data/macaque_cortex/h5ad/'

for file_path in tqdm(macaque1):
    adata = load_macaque_data(file_path)

    file_name = file_path.split('/')[-1]
    h5ad_filename = file_name[: file_name.rfind('round')] + 'h5ad'

    adata.write(h5ad_dic + h5ad_filename)

    del adata
    gc.collect()
