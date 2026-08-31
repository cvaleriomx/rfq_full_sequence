"""Entrenamiento, métricas y persistencia del modelo subrogado."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Sequence

import joblib
import numpy as np
from sklearn.compose import TransformedTargetRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler


def build_model(
    hidden_layers: Sequence[int] = (32, 32, 32, 32),
    seed: int = 2026,
    max_iter: int = 2000,
    log_targets: bool = True,
) -> TransformedTargetRegressor:
    """Crea una MLP con escalado independiente de entradas y salidas."""

    regressor = Pipeline(
        [
            ("scale_inputs", StandardScaler()),
            (
                "mlp",
                MLPRegressor(
                    hidden_layer_sizes=tuple(hidden_layers),
                    activation="tanh",
                    solver="adam",
                    learning_rate_init=1.0e-3,
                    batch_size="auto",
                    max_iter=max_iter,
                    early_stopping=True,
                    validation_fraction=0.15,
                    n_iter_no_change=80,
                    random_state=seed,
                ),
            ),
        ]
    )
    if log_targets:
        target_transformer = Pipeline(
            [
                (
                    "log",
                    FunctionTransformer(
                        np.log, inverse_func=np.exp, validate=True, check_inverse=False
                    ),
                ),
                ("scale_targets", StandardScaler()),
            ]
        )
    else:
        target_transformer = StandardScaler()
    return TransformedTargetRegressor(regressor=regressor, transformer=target_transformer)


def _per_output_metrics(
    y_true: np.ndarray, y_pred: np.ndarray, output_names: Sequence[str]
) -> dict[str, dict[str, float]]:
    result = {}
    for index, name in enumerate(output_names):
        true = y_true[:, index]
        predicted = y_pred[:, index]
        rmse = float(np.sqrt(mean_squared_error(true, predicted)))
        scale = float(np.ptp(true))
        if scale == 0.0:
            scale = max(float(np.mean(np.abs(true))), np.finfo(float).eps)
        result[name] = {
            "r2": float(r2_score(true, predicted)),
            "rmse": rmse,
            "nrmse_range": rmse / scale,
        }
    return result


def train_model(
    inputs: np.ndarray,
    outputs: np.ndarray,
    metadata: dict,
    hidden_layers: Sequence[int] = (32, 32, 32, 32),
    seed: int = 2026,
    max_iter: int = 2000,
    test_fraction: float = 0.2,
) -> tuple[dict, dict, np.ndarray, np.ndarray]:
    """Entrena y devuelve el paquete serializable, métricas y datos de paridad."""

    x_train, x_test, y_train, y_test = train_test_split(
        inputs, outputs, test_size=test_fraction, random_state=seed, shuffle=True
    )
    log_targets = bool(np.all(outputs > 0.0))
    model = build_model(hidden_layers, seed, max_iter, log_targets=log_targets)
    model.fit(x_train, y_train)
    predictions = model.predict(x_test)
    output_names = metadata["output_names"]
    per_output = _per_output_metrics(y_test, predictions, output_names)
    metrics = {
        "n_train": len(x_train),
        "n_test": len(x_test),
        "r2_global_variance_weighted": float(
            r2_score(y_test, predictions, multioutput="variance_weighted")
        ),
        "outputs": per_output,
    }
    fitted_mlp = model.regressor_.named_steps["mlp"]
    metrics["epochs"] = int(fitted_mlp.n_iter_)
    metrics["final_training_loss"] = float(fitted_mlp.loss_)
    metrics["log_targets"] = log_targets

    bundle = {
        "model": model,
        "dataset_metadata": metadata,
        "training": {
            "hidden_layers": list(hidden_layers),
            "seed": seed,
            "max_iter": max_iter,
            "test_fraction": test_fraction,
        },
        "metrics": metrics,
    }
    return bundle, metrics, y_test, predictions


def save_model_bundle(path: str | Path, bundle: dict) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, destination)
    return destination


def load_model_bundle(path: str | Path) -> dict:
    return joblib.load(path)


def save_metrics(path: str | Path, metrics: dict) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n")
    return destination


def save_parity_plot(
    path: str | Path,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    output_names: Sequence[str],
) -> Path:
    """Guarda una figura de paridad para cada observable."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    mpl_cache = destination.parent / ".mplconfig"
    mpl_cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    import matplotlib.pyplot as plt

    n_outputs = y_true.shape[1]
    n_columns = min(3, n_outputs)
    n_rows = int(np.ceil(n_outputs / n_columns))
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(4.3 * n_columns, 4.0 * n_rows))
    axes = np.atleast_1d(axes).ravel()
    for index, axis in enumerate(axes):
        if index >= n_outputs:
            axis.set_visible(False)
            continue
        true = y_true[:, index]
        predicted = y_pred[:, index]
        lower = min(float(true.min()), float(predicted.min()))
        upper = max(float(true.max()), float(predicted.max()))
        axis.scatter(true, predicted, s=14, alpha=0.75)
        axis.plot((lower, upper), (lower, upper), "k--", linewidth=1)
        axis.set_title(output_names[index])
        axis.set_xlabel("Simulación")
        axis.set_ylabel("Subrogado")
        axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(destination, dpi=160)
    plt.close(figure)
    return destination
