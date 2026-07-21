#!/usr/bin/env python3
"""Inspect ExoCTK data archive metadata without installing or extracting it."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from email.message import Message
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys
import tarfile
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urlparse
from urllib.request import Request, urlopen
import zipfile


ARCHIVE_URL_RE = re.compile(
    r"https?://[^\s\"'<>]+(?:\.tar\.gz|\.tgz|\.zip|"
    r"/shared/static/[A-Za-z0-9]+)"
)
DEFAULT_MAX_DOWNLOAD_BYTES = 100 * 1024 * 1024
CHUNK_SIZE = 1024 * 1024
IGNORED_DIRECTORY_NAMES = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".tox",
    "__pycache__",
    "build",
    "dist",
}
SAFE_HTTP_HEADERS = {
    "accept-ranges",
    "content-disposition",
    "content-length",
    "content-range",
    "content-type",
    "etag",
    "last-modified",
    "location",
}


def _patch_version(repo_root: Path) -> str | None:
    """Return the ``vMAJOR.MINOR`` data suffix used by ``utils.py``."""

    version_file = repo_root / "exoctk" / "_version.py"
    try:
        match = re.search(
            r"__version__\s*=\s*['\"]([^'\"]+)",
            version_file.read_text(encoding="utf-8"),
        )
    except OSError:
        return None
    if not match:
        return None
    parts = match.group(1).split(".")
    return "v" + ".".join(parts[:2]) if len(parts) >= 2 else None


def discover_archive_urls(repo_root: Path) -> dict[str, list[str]]:
    """Find archive URLs and the repository files that reference them."""

    references: dict[str, set[str]] = defaultdict(set)
    patch_version = _patch_version(repo_root)
    for path in repo_root.rglob("*"):
        if (
            not path.is_file()
            or IGNORED_DIRECTORY_NAMES.intersection(path.parts)
        ):
            continue
        try:
            if path.stat().st_size > 5 * 1024 * 1024:
                continue
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for match in ARCHIVE_URL_RE.finditer(text):
            url = match.group(0).rstrip(".,;)]}")
            if "{PATCHVER}" in url and patch_version:
                url = url.replace("{PATCHVER}", patch_version)
            references[url].add(str(path.relative_to(repo_root)))
    return {url: sorted(paths) for url, paths in sorted(references.items())}


def _headers_dict(headers: Message) -> dict[str, str]:
    return {
        key.lower(): value
        for key, value in headers.items()
        if key.lower() in SAFE_HTTP_HEADERS
    }


def _safe_url(url: str) -> str:
    """Remove query strings and redact Box's ephemeral signed paths."""

    parsed = urlparse(url)
    if parsed.hostname and parsed.hostname.endswith("boxcloud.com"):
        path = "/<redacted-signed-download-path>"
    else:
        path = parsed.path
    sanitized = parsed._replace(
        path=path, params="", query="", fragment=""
    )
    return sanitized.geturl()


def _total_size(headers: dict[str, str]) -> int | None:
    content_range = headers.get("content-range", "")
    match = re.search(r"/(\d+)$", content_range)
    if match:
        return int(match.group(1))
    content_length = headers.get("content-length")
    if content_length and content_length.isdigit():
        return int(content_length)
    return None


def _filename(headers: dict[str, str], final_url: str) -> str | None:
    disposition = headers.get("content-disposition", "")
    utf8_match = re.search(r"filename\*=UTF-8''([^;]+)", disposition, re.I)
    plain_match = re.search(r'filename="?([^";]+)', disposition, re.I)
    if utf8_match:
        return Path(unquote(utf8_match.group(1))).name
    if plain_match:
        return Path(plain_match.group(1)).name
    name = Path(urlparse(final_url).path).name
    return name or None


def _request_metadata(url: str, method: str, timeout: float):
    headers = {"User-Agent": "ExoCTK-data-audit/1"}
    if method == "GET":
        headers["Range"] = "bytes=0-0"
    request = Request(url, headers=headers, method=method)
    try:
        response = urlopen(request, timeout=timeout)
    except HTTPError as exc:
        return (
            exc.code,
            _safe_url(exc.geturl()),
            _headers_dict(exc.headers),
            str(exc),
        )
    with response:
        if method == "GET":
            response.read(1)
        return (
            response.status,
            _safe_url(response.geturl()),
            _headers_dict(response.headers),
            None,
        )


def inspect_http_metadata(url: str, timeout: float = 30.0) -> dict:
    """Follow redirects and report headers, with a one-byte GET fallback."""

    attempts = []
    for method in ("HEAD", "GET"):
        try:
            status, final_url, headers, error = _request_metadata(
                url, method, timeout
            )
        except URLError as exc:
            attempts.append({"method": method, "error": str(exc)})
            continue
        attempt = {
            "method": method,
            "status": status,
            "final_url": final_url,
            "headers": headers,
        }
        if error:
            attempt["error"] = error
        attempts.append(attempt)
        if 200 <= status < 400:
            return {
                "url": url,
                "status": status,
                "final_url": final_url,
                "resolved_filename": _filename(headers, final_url),
                "compressed_size_bytes": _total_size(headers),
                "headers": headers,
                "probe_method": method,
                "attempts": attempts,
            }
    detail = attempts[-1].get("error", "no successful HTTP response")
    raise RuntimeError(f"Unable to access {url}: {detail}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_archive(
    url: str,
    output_dir: Path,
    metadata: dict,
    max_download_bytes: int = DEFAULT_MAX_DOWNLOAD_BYTES,
    allow_large: bool = False,
    timeout: float = 60.0,
) -> Path:
    """Download an explicitly requested archive, but never extract it."""

    size = metadata.get("compressed_size_bytes")
    if size is not None and size > max_download_bytes and not allow_large:
        raise RuntimeError(
            f"Refusing {size}-byte download; pass --allow-large after "
            "reviewing "
            "disk, runtime, and extraction requirements"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = metadata.get("resolved_filename") or Path(
        urlparse(url).path
    ).name
    if not filename:
        raise RuntimeError(f"Could not determine a filename for {url}")
    destination = output_dir / Path(filename).name
    request = Request(url, headers={"User-Agent": "ExoCTK-data-audit/1"})
    written = 0
    with (
        urlopen(request, timeout=timeout) as response,
        destination.open("wb") as out,
    ):
        while True:
            chunk = response.read(CHUNK_SIZE)
            if not chunk:
                break
            written += len(chunk)
            if written > max_download_bytes and not allow_large:
                out.close()
                destination.unlink(missing_ok=True)
                raise RuntimeError(
                    "Download exceeded the safety limit; partial file removed"
                )
            out.write(chunk)
    return destination


def _member_report(names: Iterable[str], sizes: Iterable[int]) -> dict:
    names = list(names)
    sizes = list(sizes)
    counts = Counter(names)
    unsafe = []
    metadata_artifacts = []
    top_level = Counter()
    for name in names:
        normalized = name.replace("\\", "/")
        pure = PurePosixPath(normalized)
        parts = pure.parts
        if parts:
            top_level[parts[0]] += 1
        basename = pure.name
        if (
            "__MACOSX" in parts
            or basename == ".DS_Store"
            or basename.startswith("._")
        ):
            metadata_artifacts.append(name)
        if (
            pure.is_absolute()
            or ".." in parts
            or re.match(r"^[A-Za-z]:", normalized)
        ):
            unsafe.append(name)
    return {
        "member_count": len(names),
        "uncompressed_size_bytes": sum(sizes),
        "top_level_paths": dict(sorted(top_level.items())),
        "duplicate_member_paths": sorted(
            name for name, count in counts.items() if count > 1
        ),
        "metadata_artifacts": metadata_artifacts,
        "unsafe_member_paths": unsafe,
    }


def inspect_archive(path: Path) -> dict:
    """Inspect ZIP or tar members without extracting them."""

    path = path.resolve()
    report = {
        "path": str(path),
        "compressed_size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            members = archive.infolist()
            report.update(
                _member_report(
                    (member.filename for member in members),
                    (member.file_size for member in members),
                )
            )
            report["archive_format"] = "zip"
            report["member_compressed_size_bytes"] = sum(
                member.compress_size for member in members
            )
        return report
    try:
        with tarfile.open(path, mode="r:*") as archive:
            members = archive.getmembers()
            report.update(
                _member_report(
                    (member.name for member in members),
                    (member.size for member in members),
                )
            )
            report["archive_format"] = "tar"
        return report
    except tarfile.TarError as exc:
        raise ValueError(
            f"Unsupported or invalid archive {path}: {exc}"
        ) from exc


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--discover",
        action="store_true",
        help="discover archive URLs in the repo",
    )
    parser.add_argument(
        "--url", action="append", default=[], help="archive URL"
    )
    parser.add_argument(
        "--archive",
        action="append",
        type=Path,
        default=[],
        help="local archive",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="download URL archives and inspect",
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--max-download-bytes",
        type=int,
        default=DEFAULT_MAX_DOWNLOAD_BYTES,
    )
    parser.add_argument("--allow-large", action="store_true")
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--json", type=Path, dest="json_path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.download and args.output_dir is None:
        raise SystemExit("--download requires --output-dir")
    if not (args.discover or args.url or args.archive):
        args.discover = True

    report = {
        "schema_version": 1,
        "repository_root": str(args.repo_root.resolve()),
        "discovered_urls": {},
        "urls": [],
        "archives": [],
        "errors": [],
    }
    urls = list(args.url)
    if args.discover:
        discovered = discover_archive_urls(args.repo_root.resolve())
        report["discovered_urls"] = discovered
        urls.extend(discovered)
    urls = list(dict.fromkeys(urls))

    for url in urls:
        try:
            metadata = inspect_http_metadata(url, timeout=args.timeout)
            metadata["references"] = report["discovered_urls"].get(url, [])
            report["urls"].append(metadata)
            if args.download:
                path = download_archive(
                    url,
                    args.output_dir,
                    metadata,
                    max_download_bytes=args.max_download_bytes,
                    allow_large=args.allow_large,
                    timeout=max(args.timeout, 60.0),
                )
                archive_report = inspect_archive(path)
                archive_report["source_url"] = url
                report["archives"].append(archive_report)
        # Report all URLs rather than stopping at the first.
        except Exception as exc:
            report["errors"].append({"source": url, "error": str(exc)})

    for path in args.archive:
        try:
            report["archives"].append(inspect_archive(path))
        except Exception as exc:
            report["errors"].append({"source": str(path), "error": str(exc)})

    output = json.dumps(report, indent=2, sort_keys=True)
    if args.json_path:
        args.json_path.write_text(output + "\n", encoding="utf-8")
    print(output)
    return 1 if report["errors"] else 0


if __name__ == "__main__":
    sys.exit(main())
