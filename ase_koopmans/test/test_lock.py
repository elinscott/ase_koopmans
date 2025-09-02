def test_lock():
    """Test timeout on Lock.acquire()."""
    from ase_koopmans.utils import Lock
    from ase_koopmans.test import must_raise

    lock = Lock('lockfile', timeout=0.3)
    with lock:
        with must_raise(TimeoutError):
            with lock:
                ...
