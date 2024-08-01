for DZ in {CD,UC}; do
    for SEC in {A,B,C,D}; do
        python 1-convert-to-anndata.py \
        --spatial_path /mnt/data5/spatial/data/colon/$DZ/$SEC/outs \
        --slide_path /mnt/data5/spatial/data/colon/$DZ/$SEC/$SEC.tif \
        --tile_width 86 \
        --output_h5ad data/$DZ-$SEC-src.h5ad
        python 2-extract-expression-features.py \
        --input_h5ad data/$DZ-$SEC-src.h5ad \
        --output_h5ad data/$DZ-$SEC-expr.h5ad \
        --model uce_4 \
        --species human
        python 3-extract-histology-features.py \
        --input_h5ad data/$DZ-$SEC-src.h5ad \
        --output_h5ad data/$DZ-$SEC-hist.h5ad
        python 4-combine-data.py \
        --source_h5ad data/$DZ-$SEC-src.h5ad \
        --expr_h5ad data/$DZ-$SEC-expr.h5ad \
        --hist_h5ad data/$DZ-$SEC-hist.h5ad \
        --output_h5ad data/$DZ-$SEC.h5ad
        rm data/*-*-*.h5ad
    done
done
