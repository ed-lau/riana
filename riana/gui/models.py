# -*- coding: utf-8 -*-

"""Shared Qt item-models for the GUI tabs."""

from __future__ import annotations

import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt


class DataFrameTableModel(QAbstractTableModel):
    """A read-only Qt table model over a pandas DataFrame (results views)."""

    def __init__(self, df: pd.DataFrame | None = None) -> None:
        super().__init__()
        self._df = df if df is not None else pd.DataFrame()

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self.beginResetModel()
        self._df = df.reset_index(drop=True)
        self.endResetModel()

    @property
    def dataframe(self) -> pd.DataFrame:
        return self._df

    def rowCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._df)

    def columnCount(self, parent: QModelIndex = QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self._df.columns)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or role != Qt.ItemDataRole.DisplayRole:
            return None
        value = self._df.iat[index.row(), index.column()]
        if isinstance(value, float):
            return f"{value:.4g}"
        return str(value)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return str(self._df.columns[section])
        return str(section)
