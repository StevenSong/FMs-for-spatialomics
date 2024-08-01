import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import tab10, viridis
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from PIL import Image


def plot_overlay(
    *,
    ax,
    cmap,
    im: Image.Image,
    spot_df: pd.DataFrame,
    spot_radius: float,
    spot_values: np.ndarray,
    mark_idx: int | None = None,
    title: str | None = None,
):
    ax.imshow(im)
    cs = []
    colors = []
    for i, (x, y) in enumerate(zip(spot_df["x"], spot_df["y"])):
        spot_value = spot_values[i]
        color = cmap(spot_value)
        if i == mark_idx:
            color = "red"

        c = Circle((x, y), spot_radius)
        cs.append(c)

        colors.append(color)

    patches = PatchCollection(cs, color=colors, alpha=1)
    ax.add_collection(patches)
    ax.set_title(title)


def qclip_norm(x, qmin=0.1, qmax=0.9, eps=1e-5):
    qmin_val, qmax_val = np.quantile(x, [qmin, qmax])
    x = np.clip(x, qmin_val, qmax_val)
    x = (x - qmin_val) / (qmax_val - qmin_val + eps)
    return x


def plot_dist_rel_center(
    *,  # enforce kwargs
    our_embed_dists: np.ndarray,
    expr_embed_dists: np.ndarray,
    imags_embed_dists: np.ndarray,
    center_idx: int,
    im: Image.Image,
    spot_df: pd.DataFrame,
    spot_radius: float,
    cmap=viridis,
    quantize: int | None = None,
    even_sized: bool = True,  # or False for even_spaced
) -> None:
    # our_dists = qclip_norm(embed_dists[center_idx])
    our_dists = our_embed_dists[center_idx]  # range[0,2] (cos)
    our_dists = our_dists / 2  # range[0,1]

    expr_dists = qclip_norm(expr_embed_dists[center_idx])
    imag_dists = qclip_norm(imags_embed_dists[center_idx])

    if quantize is not None:
        if even_sized:
            bin_fn = pd.qcut
        else:  # even_spaced
            bin_fn = pd.cut
        our_clusters = bin_fn(our_dists, quantize, labels=range(quantize))
        expr_clusters = bin_fn(expr_dists, quantize, labels=range(quantize))
        imag_clusters = bin_fn(imag_dists, quantize, labels=range(quantize))

        our_dists = our_clusters
        expr_dists = expr_clusters
        imag_dists = imag_clusters

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    axs[0, 1].imshow(im)
    ref = spot_df.iloc[center_idx]
    x, y = ref["x"], ref["y"]
    axs[0, 1].add_patch(Circle((x, y), spot_radius, color="red"))

    plot_overlay(
        ax=axs[0, 0],
        im=im,
        spot_df=spot_df,
        spot_radius=spot_radius,
        spot_values=our_dists,
        mark_idx=center_idx,
        cmap=cmap,
        title="Ours",
    )
    plot_overlay(
        ax=axs[1, 0],
        im=im,
        spot_df=spot_df,
        spot_radius=spot_radius,
        spot_values=expr_dists,
        mark_idx=center_idx,
        cmap=cmap,
        title="Expression Embeddings",
    )
    plot_overlay(
        ax=axs[1, 1],
        im=im,
        spot_df=spot_df,
        spot_radius=spot_radius,
        spot_values=imag_dists,
        mark_idx=center_idx,
        cmap=cmap,
        title="Image Embeddings",
    )

    plt.show()


def plot_clusters(
    *,  # enforce kwargs
    im: Image.Image,
    spot_df: pd.DataFrame,
    spot_radius: float,
    learned_clusters: np.ndarray,
    expr_clusters: np.ndarray,
    imag_clusters: np.ndarray,
    cmap=tab10,
):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
    ax2.imshow(im)
    plot_overlay(
        ax=ax1,
        im=im,
        spot_df=spot_df,
        spot_radius=spot_radius,
        spot_values=learned_clusters,
        cmap=cmap,
        title="Ours",
    )
    plot_overlay(
        ax=ax3,
        im=im,
        spot_df=spot_df,
        spot_radius=spot_radius,
        spot_values=expr_clusters,
        cmap=cmap,
        title="Expression Embeddings",
    )
    plot_overlay(
        ax=ax4,
        im=im,
        spot_df=spot_df,
        spot_radius=spot_radius,
        spot_values=imag_clusters,
        cmap=cmap,
        title="Image Embeddings",
    )
    plt.show()
