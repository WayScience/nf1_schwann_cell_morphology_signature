import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.utils import resample
from typing import Union


# Define function to preform bootstrapping
def bootstrap_roc_auc(
    y_true: Union[np.ndarray, list],
    y_pred: Union[np.ndarray, list],
    n_bootstraps: int = 1000,
    random_seed: int = 0,
) -> np.ndarray:
    """
    Perform bootstrapping to compute the distribution of ROC AUC scores.

    Parameters:
    ----------
    y_true : array-like of shape (n_samples,)
        True binary labels (0 or 1) for the dataset.

    y_pred : array-like of shape (n_samples,)
        Predicted probabilities or scores for the positive class.

    n_bootstraps : int, optional, default=1000
        Number of bootstrap iterations to perform.

    random_seed : int, optional, default=0
        Random seed for reproducibility.

    Returns:
    -------
    bootstrapped_scores : np.ndarray
        An array of bootstrapped ROC AUC scores.
    """
    rng = np.random.default_rng(random_seed)  # Create a reproducible random generator
    bootstrapped_scores = []

    for _ in range(n_bootstraps):
        indices = rng.choice(
            len(y_true), size=len(y_true), replace=True
        )  # Reproducible resampling
        if len(np.unique(y_true[indices])) < 2:
            continue
        score = roc_auc_score(y_true[indices], y_pred[indices])
        bootstrapped_scores.append(score)

    return np.array(bootstrapped_scores)
