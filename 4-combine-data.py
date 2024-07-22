import argparse

from anndata import read_h5ad


def combine(
    *,  # enforce kwargs
    source_h5ad: str,
    expr_h5ad: str,
    hist_h5ad: str,
    output_h5ad: str,
):
    src = read_h5ad(source_h5ad)
    expr = read_h5ad(expr_h5ad)
    hist = read_h5ad(hist_h5ad)

    intersect_barcodes = src.obs.index.intersection(expr.obs.index).intersection(
        hist.obs.index
    )

    # Ensures that after filtering, obs are in the same order
    assert (
        src.obs.index.is_monotonic_increasing and not src.obs.index.duplicated().any()
    )
    assert (
        expr.obs.index.is_monotonic_increasing and not expr.obs.index.duplicated().any()
    )
    assert (
        hist.obs.index.is_monotonic_increasing and not hist.obs.index.duplicated().any()
    )

    src_mask = src.obs.index.isin(intersect_barcodes)
    expr_mask = expr.obs.index.isin(intersect_barcodes)
    hist_mask = hist.obs.index.isin(intersect_barcodes)

    filtered = src[src_mask].copy()
    filtered.obsm["X_expr"] = expr[expr_mask].obsm["X_expr"]
    filtered.obsm["X_hist"] = hist[hist_mask].obsm["X_hist"]

    filtered.write_h5ad(output_h5ad)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source_h5ad",
        required=True,
        help="Path to source spatial AnnData.",
    )
    parser.add_argument(
        "--expr_h5ad",
        required=True,
        help="Path to extracted expression AnnData.",
    )
    parser.add_argument(
        "--hist_h5ad",
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
        source_h5ad=args.source_h5ad,
        expr_h5ad=args.expr_h5ad,
        hist_h5ad=args.hist_h5ad,
        output_h5ad=args.output_h5ad,
    )
