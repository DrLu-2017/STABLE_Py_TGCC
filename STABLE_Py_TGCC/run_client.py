"""
Entry point for launching the STABLE_Py_TGCC Desktop Client GUI.
Usage:
    python run_client.py
"""

from __future__ import annotations
import sys
from PySide6.QtCore import QLocale
from PySide6.QtWidgets import QApplication, QDialog

from stable_tgcc.auth.user_vault import UserVaultManager
from stable_tgcc.client.app_icon import set_windows_app_id, get_app_icon
from stable_tgcc.client.busy_overlay import StartupSplashDialog
from stable_tgcc.client.login_dialog import LoginDialog
from stable_tgcc.client.main_window import STABLEMainWindow


def main() -> int:
    # Must be called BEFORE creating QApplication so Windows 10/11 taskbar uses our custom .ico
    set_windows_app_id("SOLEIL.STABLE.PyTGCC.A100.2.0")

    app = QApplication(sys.argv)
    # Force '.' decimal separator everywhere (avoid ',' on French/European Windows locales)
    c_loc = QLocale.c()
    c_loc.setNumberOptions(QLocale.NumberOption.OmitGroupSeparator)
    QLocale.setDefault(c_loc)
    app.setStyle("Fusion")
    icon = get_app_icon()
    app.setWindowIcon(icon)

    vault = UserVaultManager()
    dlg = LoginDialog(vault)
    dlg.setWindowIcon(icon)
    if dlg.exec() != QDialog.DialogCode.Accepted or dlg.session is None:
        return 0

    # Show animated startup splash with spinning ring & progress bar immediately after login
    splash = StartupSplashDialog(username=dlg.session.username)
    splash.show()
    app.processEvents()

    win = STABLEMainWindow(vault=vault, session=dlg.session, startup_cb=splash.set_step)
    win.setWindowIcon(icon)

    splash.set_step(100, "Workspace ready — launching main window...")
    splash.close()

    win.show()
    app.processEvents()
    win.post_show_auto_connect()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
