import sys, os, io


class captured_stdout:
    def __init__(self, filename):
        self.old_stdout = None
        self.fnull = None
        self._filename = filename

    def __enter__(self):
        self.F = open(self._filename, 'w')
        try:
            self.old_error = os.dup(sys.stderr.fileno())
            os.dup2(self.F.fileno(), sys.stderr.fileno())
        except (AttributeError, io.UnsupportedOperation):
            self.old_error = None
        return self.F

    def __exit__(self, exc_type, exc_value, traceback):
        if self.old_error is not None:
            os.dup2(self.old_error, sys.stderr.fileno())

        self.F.close()