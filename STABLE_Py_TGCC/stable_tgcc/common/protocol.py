"""
NDJSON Streaming & JSON-RPC Protocol helper between STABLE Client and TGCC Server.
"""

from __future__ import annotations
import json
from typing import Dict, Any


def encode_ndjson_event(event_type: str, payload: Dict[str, Any]) -> bytes:
    """Encode an event dictionary into a single NDJSON line (UTF-8 bytes)."""
    msg = {"event": event_type, "data": payload}
    return (json.dumps(msg, separators=(",", ":")) + "\n").encode("utf-8")


def decode_ndjson_line(line: str | bytes) -> Dict[str, Any]:
    """Decode a single NDJSON line from the server stream."""
    if isinstance(line, bytes):
        line = line.decode("utf-8", errors="replace")
    return json.loads(line.strip())
