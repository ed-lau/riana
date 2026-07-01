# -*- coding: utf-8 -*-

"""Shared Qt item-models for the GUI tabs."""

from __future__ import annotations

import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt


class DataFrameTableModel(QAbstractTableModel):
    """A read-only Qt table model over a pandas DataFrame (results views).

    ``data()`` is on the render hot path: Qt calls it once per visible cell on
    every repaint (scroll, hover, selection change). Reading each cell with
    ``DataFrame.iat`` carries pandas' scalar-lookup overhead (~7 µs/cell), so a
    single screenful repaint spends 4–13 ms just fetching text — which makes a
    large results frame (a multi-file ``integrate`` concat is 10⁵–10⁶ rows) feel
    unresponsive under scrolling and selection.

    So each column is cached as a native-dtype numpy array once per (re)assignment
    (:meth:`_rebuild_cache`) and ``data()`` indexes those arrays — ~25–80× faster
    per cell, with the rendered string byte-identical to the old ``iat`` path. The
    cache is **per column, not one 2-D object array**, so every column keeps its
    dtype (floats stay float for the ``:.4g`` format) and there is no whole-frame
    stringify — that would cost seconds on a 10⁶-row frame and defeat the point.
    Building the cache is O(columns) (each :meth:`~pandas.Series.to_numpy` is a
    cheap view), so it adds negligibly to the model reset it rides along with.
    """

    def __init__(self, df: pd.DataFrame | None = None) -> None:
        super().__init__()
        self._df = df if df is not None else pd.DataFrame()
        self._cols: list = []
        self._rebuild_cache()

    def _rebuild_cache(self) -> None:
        """Cache each column as a native-dtype ndarray for fast ``data()`` reads.

        Indexed by **position** (``iloc[:, i]``) so it stays correct even if two
        columns share a name, and kept per-column so the element dtype — and thus
        the float-vs-str formatting decision in ``data()`` — matches what
        ``DataFrame.iat`` returned exactly.
        """
        self._cols = [self._df.iloc[:, i].to_numpy()
                      for i in range(self._df.shape[1])]

    def set_dataframe(self, df: pd.DataFrame) -> None:
        self.beginResetModel()
        self._df = df.reset_index(drop=True)
        self._rebuild_cache()
        self.endResetModel()

    def sort(self, column: int, order: Qt.SortOrder = Qt.SortOrder.AscendingOrder) -> None:
        """Reorder the backing frame by *column* (called on header click).

        Sorting happens in pandas (vectorised, C-level) and mutates the model's
        own ``_df`` in place — there is no proxy. Row-click handlers therefore
        keep mapping a view row straight through ``dataframe.iloc[row]`` with no
        ``mapToSource`` translation, and the sort scales to large result frames
        where a ``QSortFilterProxyModel`` doing per-cell Python comparisons would
        not. Sorting on the raw column also keeps numeric order numeric, instead
        of the lexicographic order the ``:.4g`` display strings would give.
        """
        if self._df.empty or not (0 <= column < len(self._df.columns)):
            return
        col = self._df.columns[column]
        ascending = order == Qt.SortOrder.AscendingOrder
        self.beginResetModel()
        # mergesort is stable, so equal keys keep their prior (e.g. load) order.
        self._df = self._df.sort_values(
            by=col, ascending=ascending, kind="mergesort", na_position="last"
        ).reset_index(drop=True)
        self._rebuild_cache()
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
        # Read from the per-column ndarray cache (not ``self._df.iat``): the value
        # and its dtype are identical, but the lookup avoids pandas' per-scalar
        # overhead on this per-cell-per-repaint hot path.
        value = self._cols[index.column()][index.row()]
        if isinstance(value, float):
            return f"{value:.4g}"
        return str(value)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role != Qt.ItemDataRole.DisplayRole:
            return None
        if orientation == Qt.Orientation.Horizontal:
            return str(self._df.columns[section])
        return str(section)
