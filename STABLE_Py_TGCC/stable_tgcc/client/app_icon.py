"""
Custom Windows `.ico` Generator & Taskbar AppUserModelID Configurator for STABLE_Py_TGCC.

Why Windows shows the default Python icon on the taskbar unless configured:
  1. Windows groups any script launched via `python.exe` under the Python executable's
     AppUserModelID unless `shell32.SetCurrentProcessExplicitAppUserModelID(...)` is called.
  2. This module sets the explicit AppUserModelID AND generates a crisp multi-resolution
     Windows `.ico` file (`assets/stable_icon.ico` with 256, 64, 48, 32, 16 px layers)
     featuring a Synchrotron Storage Ring, RF Cavity Wave, and Electron Bunches.
"""

from __future__ import annotations
from pathlib import Path
import ctypes
import math
import os
import struct

from PySide6.QtCore import Qt, QByteArray, QBuffer, QIODevice, QPointF, QRectF
from PySide6.QtGui import (
    QIcon,
    QImage,
    QPainter,
    QPen,
    QBrush,
    QColor,
    QLinearGradient,
    QRadialGradient,
    QPainterPath,
    QFont,
)
from PySide6.QtWidgets import QApplication


def set_windows_app_id(app_id: str = "SOLEIL.STABLE.PyTGCC.A100.2.0") -> None:
    """Register explicit Windows AppUserModelID so the Windows 10/11 taskbar shows our custom `.ico`."""
    if os.name == "nt":
        try:
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
        except Exception:
            pass


def _render_icon_image(size: int = 256) -> QImage:
    """Render the STABLE synchrotron ring + RF wave + electron bunch icon at `size x size`."""
    img = QImage(size, size, QImage.Format.Format_ARGB32)
    img.fill(Qt.GlobalColor.transparent)

    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    p.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)

    s = float(size)
    pad = s * 0.04
    rect = QRectF(pad, pad, s - 2 * pad, s - 2 * pad)

    # 1. Deep navy-indigo rounded-square badge background
    bg_grad = QLinearGradient(0, 0, s, s)
    bg_grad.setColorAt(0.0, QColor(15, 23, 42))     # slate-900
    bg_grad.setColorAt(0.55, QColor(30, 58, 138))   # blue-900
    bg_grad.setColorAt(1.0, QColor(17, 24, 39))     # gray-900
    p.setPen(QPen(QColor(56, 189, 248, 200), max(1.5, s * 0.022)))
    p.setBrush(QBrush(bg_grad))
    p.drawRoundedRect(rect, s * 0.22, s * 0.22)

    cx, cy = s * 0.50, s * 0.46
    r_ring = s * 0.33

    # 2. Outer synchrotron storage ring glow & track
    ring_pen_glow = QPen(QColor(56, 189, 248, 75), max(2.0, s * 0.085))
    p.setPen(ring_pen_glow)
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawEllipse(QPointF(cx, cy), r_ring, r_ring)

    ring_pen = QPen(QColor(56, 189, 248), max(1.5, s * 0.032))
    p.setPen(ring_pen)
    p.drawEllipse(QPointF(cx, cy), r_ring, r_ring)

    # 3. RF Cavity Double-Harmonic Flat-Potential Waveform across center
    wave_path = QPainterPath()
    n_pts = 80
    x0 = cx - r_ring * 0.82
    x1 = cx + r_ring * 0.82
    for i in range(n_pts + 1):
        u = i / float(n_pts)
        x = x0 + u * (x1 - x0)
        phase = (u - 0.5) * 2.0 * math.pi * 1.35
        # Fundamental + 4th harmonic waveform
        y = cy - (s * 0.115) * (math.sin(phase) - 0.28 * math.sin(2.0 * phase))
        if i == 0:
            wave_path.moveTo(x, y)
        else:
            wave_path.lineTo(x, y)

    p.setPen(QPen(QColor(34, 197, 94, 235), max(1.5, s * 0.038), Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
    p.drawPath(wave_path)

    # 4. Bright Golden-Orange Electron Bunches orbiting on the ring
    for ang_deg in (35.0, 145.0, 270.0):
        ang = math.radians(ang_deg)
        bx = cx + r_ring * math.cos(ang)
        by = cy - r_ring * math.sin(ang)
        br = max(2.0, s * 0.062)
        b_grad = QRadialGradient(QPointF(bx, by), br * 1.4)
        b_grad.setColorAt(0.0, QColor(254, 240, 138))
        b_grad.setColorAt(0.55, QColor(249, 115, 22))
        b_grad.setColorAt(1.0, QColor(249, 115, 22, 0))
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QBrush(b_grad))
        p.drawEllipse(QPointF(bx, by), br * 1.4, br * 1.4)
        p.setBrush(QBrush(QColor(254, 249, 195)))
        p.drawEllipse(QPointF(bx, by), br * 0.55, br * 0.55)

    # 5. Bottom "STABLE" label banner (for sizes >= 32)
    if size >= 32:
        font = QFont("Segoe UI", max(6, int(s * 0.135)), QFont.Weight.Black)
        p.setFont(font)
        p.setPen(QColor(255, 255, 255))
        text_rect = QRectF(0, s * 0.72, s, s * 0.24)
        p.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, "STABLE")

    p.end()
    return img


def _qimage_to_png_bytes(img: QImage) -> bytes:
    ba = QByteArray()
    buf = QBuffer(ba)
    buf.open(QIODevice.OpenModeFlag.WriteOnly)
    img.save(buf, "PNG")
    return bytes(ba)


def ensure_app_ico_files() -> Tuple[Path, Path]:
    """
    Generate `assets/stable_icon.ico` (multi-size Windows ICO) and `assets/stable_icon.png`
    if not already present, and return `(ico_path, png_path)`.
    """
    assets_dir = Path(__file__).resolve().parents[2] / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)
    ico_path = assets_dir / "stable_icon.ico"
    png_path = assets_dir / "stable_icon.png"

    sizes = [256, 64, 48, 32, 16]
    png_payloads = []
    for sz in sizes:
        im = _render_icon_image(sz)
        if sz == 256:
            im.save(str(png_path), "PNG")
        png_payloads.append((sz, _qimage_to_png_bytes(im)))

    # Pack PNG layers into standard Windows ICO container format
    num_images = len(png_payloads)
    header = struct.pack("<HHH", 0, 1, num_images)
    offset = 6 + 16 * num_images
    dir_entries = []
    for sz, data in png_payloads:
        w_byte = 0 if sz >= 256 else sz
        h_byte = 0 if sz >= 256 else sz
        entry = struct.pack(
            "<BBBBHHII",
            w_byte,
            h_byte,
            0,      # palette colors
            0,      # reserved
            1,      # color planes
            32,     # bits per pixel
            len(data),
            offset,
        )
        dir_entries.append(entry)
        offset += len(data)

    with ico_path.open("wb") as f:
        f.write(header)
        for e in dir_entries:
            f.write(e)
        for _, data in png_payloads:
            f.write(data)

    return ico_path, png_path


def get_app_icon() -> QIcon:
    """Return a `QIcon` loaded from `assets/stable_icon.ico` and apply it to `QApplication`."""
    ico_path, png_path = ensure_app_ico_files()
    icon = QIcon(str(ico_path))
    icon.addFile(str(png_path))
    app = QApplication.instance()
    if isinstance(app, QApplication):
        app.setWindowIcon(icon)
    return icon
