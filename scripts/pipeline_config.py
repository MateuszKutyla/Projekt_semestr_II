#!/usr/bin/env python3
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "config.yaml"


def load_config(config_path=None):
    path = Path(config_path) if config_path else DEFAULT_CONFIG
    if not path.is_absolute():
        path = PROJECT_ROOT / path

    if not path.exists():
        return {}

    if yaml is None:
        raise RuntimeError("Do wczytania config/config.yaml wymagany jest pakiet PyYAML.")

    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def get_value(config, dotted_key, default=None):
    current = config
    for key in dotted_key.split("."):
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def project_path(value):
    path = Path(value)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path
