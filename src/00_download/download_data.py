"""
Download all files from a restricted Google Drive folder using OAuth credentials.

Usage:
    python src/00_download/download_data.py

Requirements:
    pip install google-api-python-client google-auth-httplib2 google-auth-oauthlib

Credentials:
    Place your credentials.json in config/credentials.json before running.
    On first run, a browser window will open for Google sign-in.
    A token.json will be saved to config/ for subsequent runs.
"""

import os
from pathlib import Path

from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

SCOPES = ["https://www.googleapis.com/auth/drive.readonly"]

FOLDER_ID = "1_y2Jj-R9n8_xV_ifoIbqlWGgoFYBfVg-"

ROOT = Path(__file__).resolve().parents[2]
CREDENTIALS_PATH = ROOT / "config" / "credentials.json"
TOKEN_PATH = ROOT / "config" / "token.json"
OUTPUT_DIR = ROOT / "data"


def get_credentials() -> Credentials:
    creds = None
    if TOKEN_PATH.exists():
        creds = Credentials.from_authorized_user_file(TOKEN_PATH, SCOPES)
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(CREDENTIALS_PATH, SCOPES)
            creds = flow.run_local_server(port=0)
        TOKEN_PATH.write_text(creds.to_json())
    return creds


def list_files(service, folder_id: str) -> list[dict]:
    files = []
    page_token = None
    while True:
        resp = service.files().list(
            q=f"'{folder_id}' in parents and trashed = false",
            fields="nextPageToken, files(id, name, mimeType)",
            pageToken=page_token,
        ).execute()
        files.extend(resp.get("files", []))
        page_token = resp.get("nextPageToken")
        if not page_token:
            break
    return files


def download_file(service, file_id: str, dest: Path) -> None:
    request = service.files().get_media(fileId=file_id)
    with open(dest, "wb") as f:
        downloader = MediaIoBaseDownload(f, request)
        done = False
        while not done:
            progress, done = downloader.next_chunk()
            if progress:
                print(f"  {int(progress.progress() * 100)}%", end="\r")
    print()


def main() -> None:
    assert CREDENTIALS_PATH.exists(), (
        f"credentials.json not found at {CREDENTIALS_PATH}\n"
        "See README.md → Setup: Google Drive credentials"
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    creds = get_credentials()
    service = build("drive", "v3", credentials=creds)

    print(f"Listing files in folder {FOLDER_ID} ...")
    files = list_files(service, FOLDER_ID)
    print(f"Found {len(files)} file(s)\n")

    for f in files:
        dest = OUTPUT_DIR / f["name"]
        if dest.exists():
            print(f"[skip]  {f['name']} (already exists)")
            continue
        print(f"[down]  {f['name']}")
        download_file(service, f["id"], dest)
        print(f"[done]  {f['name']}")

    print("\nAll downloads complete.")


if __name__ == "__main__":
    main()
