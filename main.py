import scanpy as sc
import sys
import matplotlib.pyplot as plt

def main():
    adata = sc.read_10x_mtx(sys.argv[1])
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    fig, ax = plt.subplots(figsize=(8,8))
    sc.pl.umap(adata, ax=ax, show=False)
    fig.savefig(sys.argv[1]+'/umap.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    main()
