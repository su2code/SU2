#!/usr/bin/env python

## \file filelock.py
#  \brief python package for filelocking
#  \author T. Lukaczyk, F. Palacios
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import os, time, errno
from random import random

# -------------------------------------------------------------------
#  File Lock Class
# -------------------------------------------------------------------
class filelock(object):
    """A file locking mechanism that has context-manager support so
    you can use it in a with statement.

    Example:
    with filelock("test.txt", timeout=2, delay=0.5):
        print("Lock acquired.")
        # Do something with the locked file

    Inputs:
        file_name - filename to lock
        timeout   - default 10sec, maximum timeout to wait for lock
        delay     - default 0.05sec, delay between each attempt to lock
                    number incremented with a random perturbation

    original source: Evan Fosmark, BSD license
    http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python/
    """

    def __init__(self, file_name, timeout=10, delay=0.05):
        """Prepare the file locker. Specify the file to lock and optionally
        the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        self.lockfile = os.path.join(os.getcwd(), "%s.lock" % file_name)
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay

    def acquire(self):
        """Acquire the lock, if possible. If the lock is in use, it check again
        every `wait` seconds. It does this until it either gets the lock or
        exceeds `timeout` number of seconds, in which case it throws
        an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT | os.O_EXCL | os.O_RDWR)
                break
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException(
                        "FileLock timeout occured for %s" % self.lockfile
                    )
                delay = self.delay * (1.0 + 0.2 * random())
                time.sleep(delay)
        self.is_locked = True

    def release(self):
        """Get rid of the lock by deleting the lockfile.
        When working in a `with` statement, this gets automatically
        called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False

    def __enter__(self):
        """Activated when used in the with statement.
        Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self

    def __exit__(self, type, value, traceback):
        """Activated at the end of the with statement.
        It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()

    def __del__(self):
        """Make sure that the FileLock instance doesn't leave a lockfile
        lying around.
        """
        self.release()


class FileLockException(Exception):
    pass


#: class filelock
