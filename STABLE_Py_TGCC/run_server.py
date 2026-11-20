"""
Entry point for starting the STABLE_Py_TGCC Remote Streaming Server on TGCC (or locally).
Usage:
    python run_server.py --host 127.0.0.1 --port 18765
"""

from __future__ import annotations
import argparse
from stable_tgcc.physics.backend import detect_available_gpus
from stable_tgcc.server.rpc_daemon import STABLERPCDaemon


def main() -> None:
    parser = argparse.ArgumentParser(description="STABLE_Py_TGCC Streaming RPC Daemon")
    parser.add_argument("--host", default="127.0.0.1", help="Bind address (default: 127.0.0.1 for SSH tunnel)")
    parser.add_argument("--port", type=int, default=18765, help="Bind port (default: 18765)")
    args = parser.parse_args()

    gpus = detect_available_gpus()
    print(f"[STABLE Server] Detected compute devices: {gpus}")
    daemon = STABLERPCDaemon(host=args.host, port=args.port)
    print(f"[STABLE Server] Listening on http://{daemon.host}:{daemon.port} ...")
    try:
        daemon.serve_forever()
    except KeyboardInterrupt:
        print("\n[STABLE Server] Shutting down.")
        daemon.shutdown()


if __name__ == "__main__":
    main()
