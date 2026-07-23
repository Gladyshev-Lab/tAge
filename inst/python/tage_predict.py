# inst/python/tage_predict.py
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional, Sequence, Union

import joblib
import pandas as pd
import sklearn.pipeline
from sklearn.exceptions import InconsistentVersionWarning
from sklearn.impute import SimpleImputer
import numpy as np

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", category=UserWarning)


PREDICTIONS_SPECIES_ADJ = {"human": 122.5, "mouse": 48, "rat": 50.4, "monkey": 39}


def _is_chronological_clock(model_path: Union[str, Path]) -> bool:
    """Return True if the model is a chronological-age clock.

    Only chronological-age clocks predict age normalised to the species
    maximum lifespan; multiplying their output by ``PREDICTIONS_SPECIES_ADJ``
    converts it back to age units (years for human, months for rodents).
    Mortality clocks output log10(hazard ratio) and normalised-age clocks
    output the normalised fraction directly, so neither must be rescaled.

    This mirrors the TACO reference application, which applies the
    species factor only to chronological clocks (``grepl("chronological")``).
    Detection is based on the model file name, following the released
    naming convention (``EN_Chronoage_...`` / ``EN_Mortality_...``).
    """
    return "chronoage" in Path(model_path).stem.lower()


def _patch_simple_imputer(imputer: SimpleImputer) -> None:
    """Patch SimpleImputer to add missing _fill_dtype attribute.

    This fixes compatibility issues when loading models trained with
    scikit-learn < 1.3.0 into environments with >= 1.3.0.

    Parameters
    ----------
    imputer : SimpleImputer
        The imputer instance to patch.
    """
    if not hasattr(imputer, "_fill_dtype"):
        if hasattr(imputer, "statistics_"):
            imputer._fill_dtype = imputer.statistics_.dtype
        else:
            imputer._fill_dtype = np.float64


def _load_clock(model_path: Union[str, Path]) -> tuple:
    """Load a clock model from file and extract feature names.

    Parameters
    ----------
    model_path : Union[str, Path]
        Path to the serialized model file.

    Returns
    -------
    tuple
        A tuple containing (clock_model, clock_genes).
    """
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError(f"Model not found: {model_path}")

    with model_path.open("rb") as f:
        clock_model = joblib.load(f)

    if isinstance(clock_model, sklearn.pipeline.Pipeline):
        for name, step in clock_model.steps:
            if isinstance(step, SimpleImputer):
                _patch_simple_imputer(step)
        clock_genes = clock_model.feature_names_in_
    else:
        clock_genes = clock_model.feature_names

    return clock_model, clock_genes


def _align_features(exprs_df: pd.DataFrame, features: Sequence[str]) -> pd.DataFrame:
    """Align expression data features with model requirements.

    Parameters
    ----------
    exprs_df : pd.DataFrame
        Expression data with features as columns.
    features : Sequence[str]
        Required feature names from the model.

    Returns
    -------
    pd.DataFrame
        Expression data with columns reordered to match model features.

    Raises
    ------
    KeyError
        If required features are missing from the expression data.
    """
    missing = [f for f in features if f not in exprs_df.columns]
    if missing:
        few = ", ".join(missing[:10])
        more = f" (+{len(missing) - 10} more)" if len(missing) > 10 else ""
        raise KeyError(f"Expression matrix is missing required features: {few}{more}")

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
    """Apply clock model to expression matrix and append predictions to annotation.

    Parameters
    ----------
    model_path : str | Path
        Path to serialized model (joblib).
    exprs_data_df : pd.DataFrame
        Normalized differential expression matrix. Columns must contain model.feature_names.
    annotation_data_df : pd.DataFrame
        Sample metadata (row count must match exprs_data_df).
    species : str
        Species name for prediction adjustment ("human", "mouse", "rat", "monkey").
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

    exprs = exprs_data_df.copy()
    exprs.columns = exprs.columns.map(str)

    cols_overlap = sum(feat in exprs.columns for feat in clock_genes)
    idx_overlap = sum(feat in exprs.index for feat in clock_genes)

    if idx_overlap > cols_overlap:
        exprs = exprs.T
        exprs.columns = exprs.columns.map(str)

    exprs.columns = exprs.columns.map(str)
    X = _align_features(exprs, clock_genes)

    if return_std:
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

    ann = annotation_data_df.copy()
    pfx = "" if prefix is None else str(prefix)
    ann.loc[:, f"{pfx}tAge"] = y
    if std is not None:
        ann.loc[:, f"{pfx}tAge_std"] = std

    # Rescale to age units only for chronological-age clocks. Mortality
    # (log10HR) and normalised-age clocks are reported on their native scale.
    if _is_chronological_clock(model_path) and species in PREDICTIONS_SPECIES_ADJ:
        ann.loc[:, f"{pfx}tAge"] = ann.loc[:, f"{pfx}tAge"] * PREDICTIONS_SPECIES_ADJ[species]

    return ann
