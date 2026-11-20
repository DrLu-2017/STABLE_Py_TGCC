"""
PySide6 User Authentication & Registration Dialog for STABLE_Py_TGCC.
Unlocks the user's encrypted TGCC SSH credentials and personal machine parameter presets.
"""

from __future__ import annotations
from typing import Optional
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QLabel,
    QLineEdit,
    QComboBox,
    QPushButton,
    QTabWidget,
    QWidget,
    QMessageBox,
)

from stable_tgcc.auth.user_vault import UserVaultManager, UserSession


class LoginDialog(QDialog):
    """
    Multi-User Login & Registration Dialog.
    Ensures every user has an isolated encrypted vault for TGCC SSH credentials and presets.
    """

    def __init__(self, vault: UserVaultManager, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.vault = vault
        self.session: Optional[UserSession] = None

        self.setWindowTitle("STABLE_Py_TGCC — User Login & Encrypted Vault")
        self.setMinimumWidth(440)
        self.setModal(True)

        root = QVBoxLayout(self)
        root.setSpacing(12)

        title = QLabel("STABLE (Python Server-Client for TGCC A100)")
        title.setFont(QFont("Segoe UI", 13, QFont.Weight.Bold))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        root.addWidget(title)

        sub = QLabel(
            "Sign in or create a local account to unlock your AES-encrypted\n"
            "TGCC SSH credentials and personal accelerator presets."
        )
        sub.setAlignment(Qt.AlignmentFlag.AlignCenter)
        sub.setStyleSheet("color: #4b5563;")
        root.addWidget(sub)

        self.tabs = QTabWidget()
        root.addWidget(self.tabs)

        # Tab 1: Login
        login_tab = QWidget()
        login_form = QFormLayout(login_tab)
        self.user_combo = QComboBox()
        self.user_combo.setEditable(True)
        existing_users = self.vault.list_usernames()
        if existing_users:
            self.user_combo.addItems(existing_users)
        self.login_pwd = QLineEdit()
        self.login_pwd.setEchoMode(QLineEdit.EchoMode.Password)
        self.login_pwd.setPlaceholderText("Enter your master password")
        self.login_pwd.returnPressed.connect(self._on_login)

        login_form.addRow("Username:", self.user_combo)
        login_form.addRow("Password:", self.login_pwd)

        btn_login = QPushButton("Unlock & Sign In")
        btn_login.setStyleSheet("background-color: #2563eb; color: white; font-weight: bold; padding: 7px;")
        btn_login.clicked.connect(self._on_login)
        login_form.addRow("", btn_login)
        self.tabs.addTab(login_tab, "Sign In")

        # Tab 2: Register New User
        reg_tab = QWidget()
        reg_form = QFormLayout(reg_tab)
        self.reg_user = QLineEdit()
        self.reg_user.setPlaceholderText("e.g., zhaolu or tgcc_login")
        self.reg_pwd1 = QLineEdit()
        self.reg_pwd1.setEchoMode(QLineEdit.EchoMode.Password)
        self.reg_pwd1.setPlaceholderText("At least 4 characters")
        self.reg_pwd2 = QLineEdit()
        self.reg_pwd2.setEchoMode(QLineEdit.EchoMode.Password)
        self.reg_pwd2.setPlaceholderText("Confirm password")
        self.reg_pwd2.returnPressed.connect(self._on_register)

        reg_form.addRow("New Username:", self.reg_user)
        reg_form.addRow("Master Password:", self.reg_pwd1)
        reg_form.addRow("Confirm Password:", self.reg_pwd2)

        btn_reg = QPushButton("Create Encrypted Account")
        btn_reg.setStyleSheet("background-color: #16a34a; color: white; font-weight: bold; padding: 7px;")
        btn_reg.clicked.connect(self._on_register)
        reg_form.addRow("", btn_reg)
        self.tabs.addTab(reg_tab, "Create Account")

        if not existing_users:
            self.tabs.setCurrentIndex(1)

    def _on_login(self) -> None:
        uname = self.user_combo.currentText().strip()
        pwd = self.login_pwd.text()
        try:
            self.session = self.vault.authenticate(uname, pwd)
            self.accept()
        except Exception as exc:
            QMessageBox.warning(self, "Authentication Failed", str(exc))

    def _on_register(self) -> None:
        uname = self.reg_user.text().strip()
        p1 = self.reg_pwd1.text()
        p2 = self.reg_pwd2.text()
        if p1 != p2:
            QMessageBox.warning(self, "Registration Error", "Passwords do not match.")
            return
        try:
            self.session = self.vault.register_user(uname, p1)
            QMessageBox.information(
                self,
                "Account Created",
                f"User '{uname}' created with default SOLEIL II & HALF presets.",
            )
            self.accept()
        except Exception as exc:
            QMessageBox.warning(self, "Registration Error", str(exc))
