import gzip
import os
import socket
import sys
import urllib.error
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import download_reads  # noqa: E402


@pytest.fixture
def no_sleep(monkeypatch):
    """Run the backoff logic without actually waiting."""
    monkeypatch.setattr(download_reads.time, "sleep", lambda _: None)


def test_dns_failure_outlasts_the_ordinary_retry_budget(no_sleep):
    """A network blip must not be mistaken for a bad accession.

    A brief outage once retired 57 good samples: every failure counted against
    the same 5-attempt budget, which a DNS outage exhausts in 25 seconds, and
    the pipeline retires whatever it cannot download.
    """
    attempts = []

    def flaky():
        attempts.append(1)
        if len(attempts) <= download_reads.MAX_RETRIES * 2:
            raise socket.gaierror(8, "nodename nor servname provided")
        return "reads"

    assert download_reads.with_retries(flaky, "flaky") == "reads"
    assert len(attempts) == download_reads.MAX_RETRIES * 2 + 1


def test_network_waiting_still_gives_up_eventually(no_sleep):
    """Waiting is bounded; an outage longer than the budget still fails."""
    def always_down():
        raise ConnectionRefusedError("connection refused")

    with pytest.raises(ConnectionRefusedError):
        download_reads.with_retries(always_down, "down")


def test_bad_accession_fails_fast(no_sleep):
    """A malformed filereport is the sample's fault: retire it promptly."""
    attempts = []

    def always_bad():
        attempts.append(1)
        raise ValueError("no fastq_ftp URLs in filereport")

    with pytest.raises(ValueError):
        download_reads.with_retries(always_bad, "bad")
    assert len(attempts) == download_reads.MAX_RETRIES


def test_a_resuming_download_is_not_charged_a_retry(no_sleep):
    """A truncated but advancing transfer is progress, not failure.

    ENA closed a 186 MB transfer early on each pass; the file climbed
    77.7 -> 77.9 -> 78.0 MB, every pass was charged against the 5-attempt
    budget, and the sample was retired while downloading correctly.
    """
    attempts = []

    def flaky_transfer():
        attempts.append(1)
        if len(attempts) <= download_reads.MAX_RETRIES * 3:
            raise download_reads.PartialTransfer(
                f"got {attempts and len(attempts) * 1000} of 186659390 bytes")
        return "complete"

    assert download_reads.with_retries(flaky_transfer, "big.fastq.gz") == "complete"
    assert len(attempts) == download_reads.MAX_RETRIES * 3 + 1


def test_a_stalled_download_is_charged_a_retry(no_sleep):
    """No forward progress is a real failure, so it must still give up."""
    attempts = []

    def stalled():
        attempts.append(1)
        raise IOError("download stalled at 0 of 186659390 bytes")

    with pytest.raises(IOError):
        download_reads.with_retries(stalled, "stalled.fastq.gz")
    assert len(attempts) == download_reads.MAX_RETRIES


def test_a_truncated_gzip_is_caught_at_download_time(tmp_path):
    """A complete-looking file whose bytes are damaged must not pass.

    One FASTQ matched its advertised Content-Length exactly and still ended
    mid-stream; nothing noticed until fastp failed and the sample was retired.
    """
    # Incompressible payload, so that lopping bytes off really does cut the
    # stream short rather than leaving a valid shorter one.
    good = tmp_path / "good.fastq.gz"
    with gzip.open(good, "wb") as fh:
        fh.write(os.urandom(256 * 1024))
    assert download_reads._gzip_ok(good)

    truncated = tmp_path / "truncated.fastq.gz"
    truncated.write_bytes(good.read_bytes()[:-1024])
    assert not download_reads._gzip_ok(truncated)

    garbage = tmp_path / "garbage.fastq.gz"
    garbage.write_bytes(b"not a gzip stream at all")
    assert not download_reads._gzip_ok(garbage)

    # An empty file decompresses to nothing without raising, so it would pass
    # a naive "did reading throw?" check.
    empty = tmp_path / "empty.fastq.gz"
    empty.write_bytes(b"")
    assert not download_reads._gzip_ok(empty)


@pytest.mark.parametrize("exc, transient", [
    (socket.gaierror(8, "nodename nor servname provided"), True),
    (urllib.error.URLError(socket.gaierror(8, "dns")), True),
    (ConnectionResetError("reset by peer"), True),
    (urllib.error.HTTPError("u", 503, "unavailable", {}, None), True),
    (urllib.error.HTTPError("u", 403, "forbidden", {}, None), True),
    (download_reads.PartialTransfer("got 1 of 2 bytes"), True),
    (urllib.error.HTTPError("u", 404, "not found", {}, None), False),
    (IOError("download stalled at 0 of 2 bytes"), False),
    (ValueError("empty filereport"), False),
])
def test_transient_and_permanent_failures_are_told_apart(exc, transient):
    assert download_reads.is_transient(exc) is transient
