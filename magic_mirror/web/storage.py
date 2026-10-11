"""SQLite-backed accounts and per-user genome storage for the web upload flow.

Deliberately dependency-free (stdlib ``sqlite3`` + ``hashlib`` + ``secrets``) so
the storage layer is testable without the optional web stack installed.  This is
a scaffold: single-file DB, bearer-token sessions, PBKDF2 password hashing —
enough to associate an uploaded genome with an account, not a hardened service.
"""

from __future__ import annotations

import hashlib
import secrets
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path

from ..sources import looks_like_23andme

_PBKDF2_ROUNDS = 200_000
_SALT_BYTES = 16
_TOKEN_BYTES = 24


class StorageError(RuntimeError):
    """A storage-level problem the API turns into a 4xx response."""


@dataclass(frozen=True, slots=True)
class GenomeRecord:
    id: int
    user_id: int
    filename: str
    stored_path: Path
    markers: int
    uploaded_at: float


def _hash_password(password: str, salt: bytes) -> bytes:
    return hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, _PBKDF2_ROUNDS)


class Storage:
    """Accounts, sessions, and genome files rooted at ``data_dir``."""

    def __init__(self, data_dir: str | Path) -> None:
        self.data_dir = Path(data_dir)
        self.genome_dir = self.data_dir / "genomes"
        self.genome_dir.mkdir(parents=True, exist_ok=True)
        self._db = sqlite3.connect(
            self.data_dir / "reflexions.db", check_same_thread=False
        )
        self._db.row_factory = sqlite3.Row
        self._db.execute("PRAGMA foreign_keys = ON")
        self._init_schema()

    def _init_schema(self) -> None:
        self._db.executescript(
            """
            CREATE TABLE IF NOT EXISTS users (
                id         INTEGER PRIMARY KEY AUTOINCREMENT,
                username   TEXT UNIQUE NOT NULL,
                pw_hash    BLOB NOT NULL,
                salt       BLOB NOT NULL,
                created_at REAL NOT NULL
            );
            CREATE TABLE IF NOT EXISTS sessions (
                token      TEXT PRIMARY KEY,
                user_id    INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                created_at REAL NOT NULL
            );
            CREATE TABLE IF NOT EXISTS genomes (
                id          INTEGER PRIMARY KEY AUTOINCREMENT,
                user_id     INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
                filename    TEXT NOT NULL,
                stored_path TEXT NOT NULL,
                markers     INTEGER NOT NULL,
                uploaded_at REAL NOT NULL
            );
            """
        )
        self._db.commit()

    def close(self) -> None:
        self._db.close()

    # --- accounts ------------------------------------------------------- #
    def create_user(self, username: str, password: str) -> int:
        username = username.strip().lower()
        if not username or not password:
            raise StorageError("username and password are required")
        salt = secrets.token_bytes(_SALT_BYTES)
        pw_hash = _hash_password(password, salt)
        try:
            cur = self._db.execute(
                "INSERT INTO users (username, pw_hash, salt, created_at) "
                "VALUES (?, ?, ?, ?)",
                (username, pw_hash, salt, time.time()),
            )
        except sqlite3.IntegrityError as exc:
            raise StorageError(f"username '{username}' is taken") from exc
        self._db.commit()
        return int(cur.lastrowid)

    def authenticate(self, username: str, password: str) -> int | None:
        row = self._db.execute(
            "SELECT id, pw_hash, salt FROM users WHERE username = ?",
            (username.strip().lower(),),
        ).fetchone()
        if row is None:
            return None
        expected = _hash_password(password, row["salt"])
        if not secrets.compare_digest(expected, row["pw_hash"]):
            return None
        return int(row["id"])

    def username_for(self, user_id: int) -> str | None:
        row = self._db.execute(
            "SELECT username FROM users WHERE id = ?", (user_id,)
        ).fetchone()
        return row["username"] if row else None

    # --- sessions ------------------------------------------------------- #
    def create_session(self, user_id: int) -> str:
        token = secrets.token_urlsafe(_TOKEN_BYTES)
        self._db.execute(
            "INSERT INTO sessions (token, user_id, created_at) VALUES (?, ?, ?)",
            (token, user_id, time.time()),
        )
        self._db.commit()
        return token

    def user_for_token(self, token: str) -> int | None:
        if not token:
            return None
        row = self._db.execute(
            "SELECT user_id FROM sessions WHERE token = ?", (token,)
        ).fetchone()
        return int(row["user_id"]) if row else None

    # --- genomes -------------------------------------------------------- #
    def save_genome(self, user_id: int, filename: str, data: bytes) -> GenomeRecord:
        """Persist an uploaded genome after a header sniff; reject non-genomes."""
        user_dir = self.genome_dir / str(user_id)
        user_dir.mkdir(parents=True, exist_ok=True)
        safe_name = Path(filename or "genome.txt").name
        stored_path = user_dir / f"{secrets.token_hex(8)}_{safe_name}"
        stored_path.write_bytes(data)

        if not looks_like_23andme(stored_path):
            stored_path.unlink(missing_ok=True)
            raise StorageError(f"{safe_name} does not look like a 23andMe raw file")

        markers = sum(
            1
            for line in data.decode("utf-8", "replace").splitlines()
            if line and not line.startswith("#") and "\t" in line
            and line.split("\t", 1)[0].strip().lower() != "rsid"
        )
        cur = self._db.execute(
            "INSERT INTO genomes (user_id, filename, stored_path, markers, uploaded_at) "
            "VALUES (?, ?, ?, ?, ?)",
            (user_id, safe_name, str(stored_path), markers, time.time()),
        )
        self._db.commit()
        return GenomeRecord(
            int(cur.lastrowid), user_id, safe_name, stored_path, markers, time.time()
        )

    def list_genomes(self, user_id: int) -> list[GenomeRecord]:
        rows = self._db.execute(
            "SELECT * FROM genomes WHERE user_id = ? ORDER BY uploaded_at DESC",
            (user_id,),
        ).fetchall()
        return [
            GenomeRecord(
                r["id"], r["user_id"], r["filename"],
                Path(r["stored_path"]), r["markers"], r["uploaded_at"],
            )
            for r in rows
        ]

    def get_genome(self, user_id: int, genome_id: int) -> GenomeRecord | None:
        """Fetch one genome, scoped to its owner so users can't read others'."""
        r = self._db.execute(
            "SELECT * FROM genomes WHERE id = ? AND user_id = ?",
            (genome_id, user_id),
        ).fetchone()
        if r is None:
            return None
        return GenomeRecord(
            r["id"], r["user_id"], r["filename"],
            Path(r["stored_path"]), r["markers"], r["uploaded_at"],
        )
