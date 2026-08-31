from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tomllib


@dataclass(frozen=True)
class PipelineConfig:
    path: Path
    root: Path
    raw: dict

    def section(self, name: str) -> dict:
        try:
            return self.raw[name]
        except KeyError as exc:
            raise KeyError(f"Falta la sección [{name}] en {self.path}") from exc

    def resolve(self, value: str | Path) -> Path:
        value = Path(value)
        return value if value.is_absolute() else self.root / value

    @property
    def output_dir(self) -> Path:
        return self.resolve(self.section("project")["output_dir"])


def load_config(path: str | Path) -> PipelineConfig:
    path = Path(path).resolve()
    with path.open("rb") as stream:
        raw = tomllib.load(stream)
    # Los archivos oficiales viven en <repo>/configs/*.toml.
    root = path.parent.parent
    return PipelineConfig(path=path, root=root, raw=raw)

