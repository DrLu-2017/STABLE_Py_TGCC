"""
Multi-User SQLite + AES/Fernet Encrypted Vault for STABLE_Py_TGCC.

Security Architecture:
  1. Each user registers with a username & master password.
  2. The master password is never stored; only PBKDF2-HMAC-SHA256 (260,000 iterations)
     salted digest is stored in SQLite (`users` table).
  3. Upon login, a per-user 256-bit Fernet key (AES-128-CBC + HMAC-SHA256) is derived in
     memory from (master_password, user_vault_salt).
  4. All sensitive TGCC SSH credentials (SSH password, private key passphrase) are encrypted
     with this per-user Fernet key before touching SQLite.
  5. Each user has isolated MachineConfig + ScanConfig presets and TGCCConnectionProfiles.
"""

from __future__ import annotations
from contextlib import closing, contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Iterator
import base64
import hashlib
import hmac
import json
import os
import sqlite3
import time

try:
    from cryptography.fernet import Fernet
    _HAS_FERNET = True
except Exception:
    Fernet = None
    _HAS_FERNET = False

from stable_tgcc.common.models import MachineConfig, ScanConfig, TGCCConnectionProfile


def _derive_fernet_key(password: str, salt: bytes) -> bytes:
    """Derive a URL-safe base64 32-byte key from password + salt using PBKDF2-HMAC-SHA256."""
    raw = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, 260_000, dklen=32)
    return base64.urlsafe_b64encode(raw)


def _encrypt_secret(secret_text: str, fernet_key: bytes) -> str:
    """Encrypt secret string with Fernet (or XOR-HMAC authenticated stream fallback)."""
    if not secret_text:
        return ""
    if _HAS_FERNET and Fernet is not None:
        f = Fernet(fernet_key)
        return f.encrypt(secret_text.encode("utf-8")).decode("ascii")
    # Fallback authenticated stream cipher using SHA256 counter mode + HMAC
    raw_key = base64.urlsafe_b64decode(fernet_key)
    nonce = os.urandom(16)
    pt = secret_text.encode("utf-8")
    ct = bytearray(len(pt))
    for i, b in enumerate(pt):
        blk = hashlib.sha256(raw_key + nonce + i.to_bytes(4, "big")).digest()
        ct[i] = b ^ blk[i % 32]
    mac = hmac.new(raw_key, nonce + bytes(ct), hashlib.sha256).digest()
    return "FALLBACK:" + base64.b64encode(nonce + mac + bytes(ct)).decode("ascii")


def _decrypt_secret(token_str: str, fernet_key: bytes) -> str:
    """Decrypt secret string encrypted by `_encrypt_secret`."""
    if not token_str:
        return ""
    if token_str.startswith("FALLBACK:"):
        raw_key = base64.urlsafe_b64decode(fernet_key)
        blob = base64.b64decode(token_str[len("FALLBACK:"):])
        nonce, mac, ct = blob[:16], blob[16:48], blob[48:]
        exp_mac = hmac.new(raw_key, nonce + ct, hashlib.sha256).digest()
        if not hmac.compare_digest(mac, exp_mac):
            raise ValueError("Invalid encryption key or tampered secret")
        pt = bytearray(len(ct))
        for i, b in enumerate(ct):
            blk = hashlib.sha256(raw_key + nonce + i.to_bytes(4, "big")).digest()
            pt[i] = b ^ blk[i % 32]
        return bytes(pt).decode("utf-8")
    if _HAS_FERNET and Fernet is not None:
        f = Fernet(fernet_key)
        return f.decrypt(token_str.encode("ascii")).decode("utf-8")
    raise RuntimeError("cryptography package required to decrypt Fernet token")


@dataclass
class UserSession:
    """Authenticated user session holding the in-memory decryption key."""
    user_id: int
    username: str
    fernet_key: bytes


class UserVaultManager:
    """
    Manages multi-user registration, authentication, encrypted TGCC SSH credentials,
    and per-user machine/scan parameter presets in SQLite.
    """

    def __init__(self, db_path: Optional[str | Path] = None):
        if db_path is None:
            base_dir = Path(__file__).resolve().parents[2] / "data"
            base_dir.mkdir(parents=True, exist_ok=True)
            db_path = base_dir / "stable_users.db"
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_schema()

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        conn = sqlite3.connect(str(self.db_path))
        conn.row_factory = sqlite3.Row
        try:
            with conn:
                yield conn
        finally:
            conn.close()

    def _init_schema(self) -> None:
        with self._connect() as conn:
            conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS users (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    username TEXT UNIQUE NOT NULL,
                    pwd_salt_hex TEXT NOT NULL,
                    pwd_hash_hex TEXT NOT NULL,
                    vault_salt_hex TEXT NOT NULL,
                    created_at REAL NOT NULL
                );

                CREATE TABLE IF NOT EXISTS user_tgcc_profiles (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id INTEGER NOT NULL,
                    profile_name TEXT NOT NULL,
                    profile_json TEXT NOT NULL,
                    encrypted_ssh_secret TEXT NOT NULL,
                    updated_at REAL NOT NULL,
                    UNIQUE(user_id, profile_name),
                    FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
                );

                CREATE TABLE IF NOT EXISTS user_presets (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    user_id INTEGER NOT NULL,
                    preset_name TEXT NOT NULL,
                    machine_json TEXT NOT NULL,
                    scan_json TEXT NOT NULL,
                    updated_at REAL NOT NULL,
                    UNIQUE(user_id, preset_name),
                    FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
                );

                CREATE TABLE IF NOT EXISTS user_last_session (
                    user_id INTEGER PRIMARY KEY,
                    machine_json TEXT NOT NULL,
                    scan_json TEXT NOT NULL,
                    updated_at REAL NOT NULL,
                    FOREIGN KEY(user_id) REFERENCES users(id) ON DELETE CASCADE
                );
                """
            )

    def list_usernames(self) -> List[str]:
        with self._connect() as conn:
            rows = conn.execute("SELECT username FROM users ORDER BY username ASC").fetchall()
            return [str(r["username"]) for r in rows]

    def register_user(self, username: str, password: str) -> UserSession:
        """Register a new user and seed default SOLEIL II & HALF presets and TGCC profile."""
        uname = username.strip()
        if not uname:
            raise ValueError("Username cannot be empty.")
        if len(password) < 4:
            raise ValueError("Password must be at least 4 characters.")

        pwd_salt = os.urandom(16)
        vault_salt = os.urandom(16)
        pwd_hash = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), pwd_salt, 260_000)

        with self._connect() as conn:
            try:
                cur = conn.execute(
                    """
                    INSERT INTO users (username, pwd_salt_hex, pwd_hash_hex, vault_salt_hex, created_at)
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    (uname, pwd_salt.hex(), pwd_hash.hex(), vault_salt.hex(), time.time()),
                )
                uid = int(cur.lastrowid)
            except sqlite3.IntegrityError:
                raise ValueError(f"User '{uname}' already exists.")

        session = UserSession(
            user_id=uid,
            username=uname,
            fernet_key=_derive_fernet_key(password, vault_salt),
        )

        # Seed default presets for this user
        self.save_user_preset(
            session,
            "SOLEIL II (Default)",
            MachineConfig.soleil2_default(),
            ScanConfig(),
        )
        self.save_user_preset(
            session,
            "HALF (Default)",
            MachineConfig.half_default(),
            ScanConfig(Df_min_kHz=20.0, Df_max_kHz=80.0, Df_step_kHz=5.0),
        )
        default_tgcc = TGCCConnectionProfile(
            profile_name="TGCC Irene A100",
            mode="ssh_tgcc",
            ssh_host="topaze.ccc.cea.fr",
            ssh_username=uname or "zhaolu",
            ccc_project="soldai",
            ccc_partition="a100",
            remote_work_dir="/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC",
            python_env_cmd="python3",
        )
        self.save_tgcc_profile(session, default_tgcc, ssh_secret="")
        return session

    def authenticate(self, username: str, password: str) -> UserSession:
        """Verify username & password and return an active `UserSession`."""
        uname = username.strip()
        with self._connect() as conn:
            row = conn.execute(
                "SELECT id, username, pwd_salt_hex, pwd_hash_hex, vault_salt_hex FROM users WHERE username = ?",
                (uname,),
            ).fetchone()
            if row is None:
                raise ValueError(f"User '{uname}' does not exist.")

            pwd_salt = bytes.fromhex(row["pwd_salt_hex"])
            expected_hash = bytes.fromhex(row["pwd_hash_hex"])
            calc_hash = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), pwd_salt, 260_000)
            if not hmac.compare_digest(calc_hash, expected_hash):
                raise ValueError("Invalid password.")

            vault_salt = bytes.fromhex(row["vault_salt_hex"])
            return UserSession(
                user_id=int(row["id"]),
                username=str(row["username"]),
                fernet_key=_derive_fernet_key(password, vault_salt),
            )

    def change_password(self, session: UserSession, old_password: str, new_password: str) -> UserSession:
        """Re-authenticate, re-encrypt all stored SSH secrets under the new password key, and update user."""
        verified = self.authenticate(session.username, old_password)
        if len(new_password) < 4:
            raise ValueError("New password must be at least 4 characters.")

        new_pwd_salt = os.urandom(16)
        new_vault_salt = os.urandom(16)
        new_pwd_hash = hashlib.pbkdf2_hmac("sha256", new_password.encode("utf-8"), new_pwd_salt, 260_000)
        new_fernet_key = _derive_fernet_key(new_password, new_vault_salt)

        with self._connect() as conn:
            rows = conn.execute(
                "SELECT id, encrypted_ssh_secret FROM user_tgcc_profiles WHERE user_id = ?",
                (verified.user_id,),
            ).fetchall()
            for r in rows:
                plain_secret = _decrypt_secret(str(r["encrypted_ssh_secret"]), verified.fernet_key)
                re_enc = _encrypt_secret(plain_secret, new_fernet_key)
                conn.execute(
                    "UPDATE user_tgcc_profiles SET encrypted_ssh_secret = ?, updated_at = ? WHERE id = ?",
                    (re_enc, time.time(), int(r["id"])),
                )
            conn.execute(
                """
                UPDATE users
                SET pwd_salt_hex = ?, pwd_hash_hex = ?, vault_salt_hex = ?
                WHERE id = ?
                """,
                (new_pwd_salt.hex(), new_pwd_hash.hex(), new_vault_salt.hex(), verified.user_id),
            )

        return UserSession(user_id=verified.user_id, username=verified.username, fernet_key=new_fernet_key)

    def save_tgcc_profile(
        self,
        session: UserSession,
        profile: TGCCConnectionProfile,
        ssh_secret: str = "",
    ) -> None:
        """Save or update a user's TGCC SSH connection profile with AES/Fernet encrypted password."""
        enc_secret = _encrypt_secret(ssh_secret, session.fernet_key)
        p_json = json.dumps(profile.to_dict(), indent=2)
        with self._connect() as conn:
            conn.execute(
                """
                INSERT INTO user_tgcc_profiles (user_id, profile_name, profile_json, encrypted_ssh_secret, updated_at)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(user_id, profile_name) DO UPDATE SET
                    profile_json = excluded.profile_json,
                    encrypted_ssh_secret = excluded.encrypted_ssh_secret,
                    updated_at = excluded.updated_at
                """,
                (session.user_id, profile.profile_name, p_json, enc_secret, time.time()),
            )

    def load_tgcc_profile(
        self,
        session: UserSession,
        profile_name: str = "TGCC Irene A100",
    ) -> Tuple[TGCCConnectionProfile, str]:
        """Load and decrypt a user's TGCC SSH connection profile and secret."""
        with self._connect() as conn:
            row = conn.execute(
                """
                SELECT profile_json, encrypted_ssh_secret
                FROM user_tgcc_profiles
                WHERE user_id = ? AND profile_name = ?
                """,
                (session.user_id, profile_name),
            ).fetchone()
            if row is None:
                prof = TGCCConnectionProfile(profile_name=profile_name, ssh_username=session.username)
                return prof, ""
            prof = TGCCConnectionProfile.from_dict(json.loads(row["profile_json"]))
            secret = _decrypt_secret(str(row["encrypted_ssh_secret"]), session.fernet_key)
            return prof, secret

    def list_user_presets(self, session: UserSession) -> List[str]:
        with self._connect() as conn:
            rows = conn.execute(
                """
                SELECT preset_name FROM user_presets
                WHERE user_id = ?
                ORDER BY CASE WHEN preset_name LIKE 'SOLEIL%' THEN 0 ELSE 1 END ASC, updated_at DESC
                """,
                (session.user_id,),
            ).fetchall()
            return [str(r["preset_name"]) for r in rows]

    def save_user_preset(
        self,
        session: UserSession,
        preset_name: str,
        machine_cfg: MachineConfig,
        scan_cfg: ScanConfig,
    ) -> None:
        """Save or update a user's named MachineConfig + ScanConfig preset."""
        machine_cfg.preset_name = preset_name
        with self._connect() as conn:
            conn.execute(
                """
                INSERT INTO user_presets (user_id, preset_name, machine_json, scan_json, updated_at)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(user_id, preset_name) DO UPDATE SET
                    machine_json = excluded.machine_json,
                    scan_json = excluded.scan_json,
                    updated_at = excluded.updated_at
                """,
                (
                    session.user_id,
                    preset_name,
                    json.dumps(machine_cfg.to_dict(), indent=2),
                    json.dumps(scan_cfg.to_dict(), indent=2),
                    time.time(),
                ),
            )

    def load_user_preset(
        self,
        session: UserSession,
        preset_name: str,
    ) -> Tuple[MachineConfig, ScanConfig]:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT machine_json, scan_json FROM user_presets WHERE user_id = ? AND preset_name = ?",
                (session.user_id, preset_name),
            ).fetchone()
            if row is None:
                raise ValueError(f"Preset '{preset_name}' not found for user '{session.username}'.")
            return (
                MachineConfig.from_dict(json.loads(row["machine_json"])),
                ScanConfig.from_dict(json.loads(row["scan_json"])),
            )

    def delete_user_preset(self, session: UserSession, preset_name: str) -> None:
        with self._connect() as conn:
            conn.execute(
                "DELETE FROM user_presets WHERE user_id = ? AND preset_name = ?",
                (session.user_id, preset_name),
            )

    def save_last_session_state(
        self,
        session: UserSession,
        machine_cfg: MachineConfig,
        scan_cfg: ScanConfig,
    ) -> None:
        """Persist the active UI MachineConfig + ScanConfig on exit so the next launch restores it."""
        with self._connect() as conn:
            conn.execute(
                """
                INSERT INTO user_last_session (user_id, machine_json, scan_json, updated_at)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(user_id) DO UPDATE SET
                    machine_json = excluded.machine_json,
                    scan_json = excluded.scan_json,
                    updated_at = excluded.updated_at
                """,
                (
                    session.user_id,
                    json.dumps(machine_cfg.to_dict(), indent=2),
                    json.dumps(scan_cfg.to_dict(), indent=2),
                    time.time(),
                ),
            )

    def load_last_session_state(
        self,
        session: UserSession,
    ) -> Optional[Tuple[MachineConfig, ScanConfig]]:
        """Load the last saved UI state from the previous session, if any."""
        with self._connect() as conn:
            row = conn.execute(
                "SELECT machine_json, scan_json FROM user_last_session WHERE user_id = ?",
                (session.user_id,),
            ).fetchone()
            if row is None:
                return None
            return (
                MachineConfig.from_dict(json.loads(row["machine_json"])),
                ScanConfig.from_dict(json.loads(row["scan_json"])),
            )
