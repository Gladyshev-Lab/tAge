# inst/python/tage_predict.py
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Union
import joblib
import pandas as pd
import sklearn.pipeline


PREDICTIONS_SPECIES_ADJ = {"human": 122.5, "mouse": 48, "rat": 50.4, "monkey": 39}


def _load_clock(model_path: Union[str, Path]):
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError(f"Model not found: {model_path}")
    with model_path.open("rb") as f:
        clock_model = joblib.load(f)

    if isinstance(clock_model, sklearn.pipeline.Pipeline):
        clock_genes = clock_model.feature_names_in_
    else:
        clock_genes = clock_model.feature_names
    
    return clock_model, clock_genes


def _align_features(exprs_df: pd.DataFrame, features: Sequence[str]) -> pd.DataFrame:
    missing = [f for f in features if f not in exprs_df.columns]
    if missing:
        few = ", ".join(missing[:10])
        more = f" (+{len(missing) - 10} more)" if len(missing) > 10 else ""
        raise KeyError(f"Expression matrix is missing required features: {few}{more}")

    # Save the order of features as in the model
    return exprs_df.loc[:, list(features)]


def predict_tAge(
    model_path: Union[str, Path],
    exprs_data_df: pd.DataFrame,
    annotation_data_df: pd.DataFrame,
    species: str,
    *,
    return_std: bool = False,
    prefix: Optional[str] = None,
) -> pd.DataFrame:
    """
    Apply clock model to expression matrix and append predictions to annotation.

    Parameters
    ----------
    model_path : str | Path
        Path to serialized model (joblib).
    exprs_data_df : pd.DataFrame
        Normalized differential expression matrix. Columns must contain model.feature_names.
    annotation_data_df : pd.DataFrame
        Sample metadata (row count must match exprs_data_df).
    return_std : bool, default False
        If True, also return predictive std (model must support return_std=True).
    prefix : str | None
        Prefix for column names (e.g., "BR_"). If None, use no prefix.

    Returns
    -------
    pd.DataFrame
        A copy of annotation_data_df with added columns:
        - f"{prefix}tAge" (always)
        - f"{prefix}tAge_std" (if return_std)
    """
    clock, clock_genes = _load_clock(model_path)

    # Copy to avoid modifying input data
    exprs = exprs_data_df.copy()
    exprs.columns = exprs.columns.map(str)

    # Auto-detect: where is more overlap of model features - in columns or in index?
    cols_overlap = sum(feat in exprs.columns for feat in clock_genes)
    idx_overlap = sum(feat in exprs.index for feat in clock_genes)

    if idx_overlap > cols_overlap:
        # It seems that genes are in the indices (rows) -> transpose
        exprs = exprs.T
        exprs.columns = exprs.columns.map(str)

    # Make sure features are strings
    exprs.columns = exprs.columns.map(str)

    X = _align_features(exprs, clock_genes)

    if return_std:
        # Check if the model supports return_std
        try:
            y, std = clock.predict(X, return_std=True)
        except TypeError as e:
            raise TypeError(
                "Model does not support return_std. Use return_std=False or "
                "supply a model with probabilistic predictions."
            ) from e
    else:
        y = clock.predict(X)
        std = None

    # Accurate copy of metadata
    ann = annotation_data_df.copy()
    pfx = "" if prefix is None else str(prefix)
    ann.loc[:, f"{pfx}tAge"] = y
    if std is not None:
        ann.loc[:, f"{pfx}tAge_std"] = std
    
    # Adjust predictions for species if applicable
    if species in PREDICTIONS_SPECIES_ADJ:
        ann.loc[:, f"{pfx}tAge"] = ann.loc[:, f"{pfx}tAge"] * PREDICTIONS_SPECIES_ADJ[species]

    return ann
