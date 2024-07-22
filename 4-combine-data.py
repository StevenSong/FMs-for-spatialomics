import argparse

from anndata import read_h5ad


def combine(
    *,  # enforce kwargs
    src_h5ad: str,
    uce_h5ad: str,
    uni_h5ad: str,
    output_h5ad: str,
):
    src = read_h5ad(src_h5ad)
    uce = read_h5ad(uce_h5ad)
    uni = read_h5ad(uni_h5ad)

    intersect_barcodes = src.obs.index.intersection(uce.obs.index).intersection(
        uni.obs.index
    )

    # Ensures that after filtering, obs are in the same order
    assert (
        src.obs.index.is_monotonic_increasing and not src.obs.index.duplicated().any()
    )
    assert (
        uce.obs.index.is_monotonic_increasing and not uce.obs.index.duplicated().any()
    )
    assert (
        uni.obs.index.is_monotonic_increasing and not uni.obs.index.duplicated().any()
    )

    filtered = src[src.obs.index.isin(intersect_barcodes)].copy()
    filtered.obsm["X_uce"] = uce[uce.obs.index.isin(intersect_barcodes)].obsm["X_uce"]
    filtered.obsm["X_uni"] = uni[uni.obs.index.isin(intersect_barcodes)].obsm["X_uni"]

    filtered.write_h5ad(output_h5ad)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--src_h5ad",
        required=True,
        help="Path to source spatial AnnData.",
    )
    parser.add_argument(
        "--uce_h5ad",
        required=True,
        help="Path to extracted expression AnnData.",
    )
    parser.add_argument(
        "--uni_h5ad",
        required=True,
        help="Path to extracted histology AnnData.",
    )
    parser.add_argument(
        "--output_h5ad",
        required=True,
        help="Path to save combined data.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    combine(
        src_h5ad=args.src_h5ad,
        uce_h5ad=args.uce_h5ad,
        uni_h5ad=args.uni_h5ad,
        output_h5ad=args.output_h5ad,
    )
