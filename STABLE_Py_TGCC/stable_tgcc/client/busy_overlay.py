"""
Animated Spinning Ring Indicator, Post-Login Startup Splash Dialog,
and Non-Blocking Async Task Progress Dialog for STABLE_Py_TGCC.
"""

from __future__ import annotations
from typing import Optional
from PySide6.QtCore import Qt, QTimer, QRectF
from PySide6.QtGui import QPainter, QColor, QPen, QFont
from PySide6.QtWidgets import (
    QDialog,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QApplication,
)

from stable_tgcc.client.app_icon import get_app_icon


class SpinningRingWidget(QWidget):
    """Smooth 60-FPS dual-arc rotating spinner widget."""

    def __init__(self, size: int = 44, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._size = size
        self.setFixedSize(size, size)
        self._angle = 0
        self._timer = QTimer(self)
        self._timer.setInterval(25)  # ~40 FPS smooth rotation
        self._timer.timeout.connect(self._tick)
        self._timer.start()

    def _tick(self) -> None:
        self._angle = (self._angle + 12) % 360
        self.update()

    def stop(self) -> None:
        self._timer.stop()

    def start(self) -> None:
        if not self._timer.isActive():
            self._timer.start()

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

        pad = 5.0
        rect = QRectF(pad, pad, self._size - 2 * pad, self._size - 2 * pad)

        # Subtle background track ring
        track_pen = QPen(QColor(226, 232, 240), 4.0)
        track_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(track_pen)
        painter.drawEllipse(rect)

        # Primary rotating arc (deep blue)
        arc_pen1 = QPen(QColor(37, 99, 235), 4.2)
        arc_pen1.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(arc_pen1)
        painter.drawArc(rect, int(-self._angle * 16), int(110 * 16))

        # Secondary counter-arc (cyan accent)
        arc_pen2 = QPen(QColor(14, 165, 233), 3.6)
        arc_pen2.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(arc_pen2)
        painter.drawArc(rect, int((-self._angle + 180) * 16), int(70 * 16))
        painter.end()


class StartupSplashDialog(QDialog):
    """
    Post-login startup progress splash window with animated spinning ring
    and step-by-step progress bar while STABLEMainWindow initializes.
    """

    def __init__(self, username: str, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setWindowTitle("Launching STABLE_Py_TGCC...")
        self.setWindowIcon(get_app_icon())
        self.setWindowFlags(
            Qt.WindowType.Dialog
            | Qt.WindowType.CustomizeWindowHint
            | Qt.WindowType.WindowTitleHint
        )
        self.setModal(True)
        self.setFixedSize(460, 175)
        self.setStyleSheet(
            "QDialog { background-color: #f8fafc; border: 1px solid #cbd5e1; }"
        )

        root = QVBoxLayout(self)
        root.setContentsMargins(22, 20, 22, 18)
        root.setSpacing(12)

        top_row = QHBoxLayout()
        top_row.setSpacing(16)
        self.spinner = SpinningRingWidget(size=48, parent=self)
        top_row.addWidget(self.spinner)

        title_col = QVBoxLayout()
        title_col.setSpacing(4)
        lbl_title = QLabel(f"Welcome, {username} — Initializing Workspace")
        lbl_title.setFont(QFont("Segoe UI", 11, QFont.Weight.Bold))
        lbl_title.setStyleSheet("color: #0f172a;")
        title_col.addWidget(lbl_title)

        self.lbl_step = QLabel("Preparing STABLE_Py_TGCC Physics & Multi-A100 Suite...")
        self.lbl_step.setFont(QFont("Segoe UI", 9))
        self.lbl_step.setStyleSheet("color: #475569;")
        title_col.addWidget(self.lbl_step)
        top_row.addLayout(title_col, stretch=1)
        root.addLayout(top_row)

        self.progress = QProgressBar(self)
        self.progress.setRange(0, 100)
        self.progress.setValue(8)
        self.progress.setTextVisible(True)
        self.progress.setFixedHeight(20)
        self.progress.setStyleSheet(
            """
            QProgressBar {
                border: 1px solid #cbd5e1;
                border-radius: 5px;
                background: #e2e8f0;
                text-align: center;
                font-weight: bold;
                color: #1e293b;
            }
            QProgressBar::chunk {
                background-color: #2563eb;
                border-radius: 4px;
            }
            """
        )
        root.addWidget(self.progress)

    def set_step(self, percent: int, message: str) -> None:
        self.progress.setValue(max(0, min(100, percent)))
        self.lbl_step.setText(message)
        self.spinner._tick()
        QApplication.processEvents()

    def closeEvent(self, event) -> None:  # noqa: N802
        self.spinner.stop()
        super().closeEvent(event)


class AsyncTaskProgressDialog(QDialog):
    """
    Non-blocking modal progress dialog with a live 60-FPS spinning ring and
    step progress bar used during SSH connection, code deployment, and batch submission.
    """

    def __init__(self, title: str, initial_text: str, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setWindowIcon(get_app_icon())
        self.setWindowFlags(
            Qt.WindowType.Dialog
            | Qt.WindowType.CustomizeWindowHint
            | Qt.WindowType.WindowTitleHint
        )
        self.setModal(True)
        self.setFixedSize(480, 165)
        self.setStyleSheet("QDialog { background-color: #ffffff; }")

        root = QVBoxLayout(self)
        root.setContentsMargins(20, 18, 20, 18)
        root.setSpacing(12)

        row = QHBoxLayout()
        row.setSpacing(16)
        self.spinner = SpinningRingWidget(size=46, parent=self)
        row.addWidget(self.spinner)

        col = QVBoxLayout()
        col.setSpacing(4)
        self.lbl_header = QLabel(title)
        self.lbl_header.setFont(QFont("Segoe UI", 10, QFont.Weight.Bold))
        self.lbl_header.setStyleSheet("color: #0f172a;")
        col.addWidget(self.lbl_header)

        self.lbl_status = QLabel(initial_text)
        self.lbl_status.setFont(QFont("Segoe UI", 9))
        self.lbl_status.setStyleSheet("color: #334155;")
        self.lbl_status.setWordWrap(True)
        col.addWidget(self.lbl_status)
        row.addLayout(col, stretch=1)
        root.addLayout(row)

        self.progress = QProgressBar(self)
        # Indeterminate marquee animation by default (range 0, 0)
        self.progress.setRange(0, 0)
        self.progress.setFixedHeight(16)
        self.progress.setStyleSheet(
            """
            QProgressBar {
                border: 1px solid #cbd5e1;
                border-radius: 4px;
                background: #f1f5f9;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #0284c7;
                border-radius: 3px;
            }
            """
        )
        root.addWidget(self.progress)

    def update_status(self, message: str, percent: Optional[int] = None) -> None:
        self.lbl_status.setText(message)
        if percent is not None:
            if self.progress.maximum() == 0:
                self.progress.setRange(0, 100)
            self.progress.setValue(max(0, min(100, percent)))

    def finish_and_close(self) -> None:
        self.spinner.stop()
        self.accept()
