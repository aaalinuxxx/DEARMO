"""Фазы docking-пайплайна ORC."""
from __future__ import annotations

INGEST = "ingest"
CLUSTER = "cluster"
PARALLEL_SEARCH = "parallel_search"
PROPERTY_SUMMARY = "property_summary"
FINALIZE = "finalize"

PHASES = (INGEST, CLUSTER, PARALLEL_SEARCH, PROPERTY_SUMMARY, FINALIZE)
