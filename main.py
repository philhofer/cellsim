import scanpy as sc
import sys
import matplotlib.pyplot as plt

N_PCS=50
N_NEIGHBORS=5

def main():
    print("load...")
    adata = sc.read_10x_mtx(sys.argv[1])
    print("lognormalize...")
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    print("pca...")
    sc.tl.pca(adata)
    print("neighbors...")
    sc.pp.neighbors(adata)
    print("umap...")
    sc.tl.umap(adata)

    print("plotting umap...")
    fig, ax = plt.subplots(figsize=(8,8))
    sc.pl.umap(adata, ax=ax, show=False)
    fig.savefig(sys.argv[1]+'/umap.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("done.")


if __name__ == "__main__":
    main()
