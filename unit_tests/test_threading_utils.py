import time

import pytest
from datetime import datetime, timedelta
import threading

from scm.plams.core.threading_utils import LimitedSemaphore


class TestLimitedSemaphore:

    def test_can_acquire_and_release_without_blocking_within_limit(self):
        # Given semaphore with max
        semaphore = LimitedSemaphore(3)

        # When acquire then release within limit
        # Then does not block
        for i in range(10):
            semaphore.acquire()
            with semaphore:
                semaphore.acquire()
                semaphore.release()
            semaphore.acquire()
            semaphore.acquire()
            semaphore.release(3)

    def test_release_above_limit_errors(self):
        # Given semaphore with max
        semaphore = LimitedSemaphore(3)

        # When release above limit
        # Then raises error
        with pytest.raises(ValueError):
            semaphore.release()
        semaphore.acquire()
        with pytest.raises(ValueError):
            semaphore.release(2)

    def test_acquire_blocks_until_released(self):
        # Given semaphore with max
        semaphore = LimitedSemaphore(3)

        # When acquire
        log_values = {}
        thread1 = self.start_in_thread(self.acquire_and_log, 0, semaphore, 3, "a1", log_values)
        thread2 = self.start_in_thread(self.acquire_and_log, 0.05, semaphore, 3, "a2", log_values)

        # Then extra threads wait for release
        time.sleep(0.1)
        semaphore.release(3)
        thread1.join()
        thread2.join()

        assert log_values["a2"] - log_values["a1"] > timedelta(seconds=0.1)

    def test_set_max_value(self):
        # Given semaphore with initial max
        semaphore = LimitedSemaphore(3)
        assert semaphore.max_value == 3

        # When increase max on thread with delay
        log_values = {}
        self.acquire_and_log(semaphore, 3, "a1", log_values)
        self.start_in_thread(self.set_max_value_and_log, 0.2, semaphore, 10, "m10", log_values)

        # Then can acquire after max increased
        self.acquire_and_log(semaphore, 7, "a2", log_values)
        assert log_values["a2"] - log_values["a1"] > timedelta(seconds=0.2)
        assert semaphore.max_value == 10

        # When decrease max
        t = self.start_in_thread(self.set_max_value_and_log, 0.0, semaphore, 2, "m2", log_values)

        # Then waits for enough releases to take effect
        self.release_and_log(semaphore, 5, "r1", log_values)
        assert semaphore.max_value == 10
        self.release_and_log(semaphore, 3, "r2", log_values)
        t.join()
        assert semaphore.max_value == 2

        # When decrease to zero
        t = self.start_in_thread(self.set_max_value_and_log, 0.0, semaphore, 0, "m0", log_values)
        self.release_and_log(semaphore, 2, "r2", log_values)
        t.join()
        assert semaphore.max_value == 0

        # Then cannot acquire until increased
        t = self.start_in_thread(self.acquire_and_log, 0, semaphore, 3, "a3", log_values)
        self.start_in_thread(self.set_max_value_and_log, 0.2, semaphore, 100, "m100", log_values)
        t.join()
        assert log_values["a3"] - log_values["m0"] > timedelta(seconds=0.2)

    def test_many_threads(self):
        # Set limit many times and perform many acquires/releases, check completes
        limits = [2, 8, 1, 16, 0, 128, 32, 0, 64, 256, 4, 0, 256, 32, 16, 0, 64]
        semaphore = LimitedSemaphore(0)
        log_values = {}

        def acquire_delay_release(n, d):
            def acquire_delay_release():
                with semaphore:
                    time.sleep(d)

            threads = [threading.Thread(target=acquire_delay_release) for _ in range(n)]
            for thread in threads:
                thread.start()
            return threads

        threads = []
        for i, limit in enumerate(limits):
            t = acquire_delay_release(128, 0.05)
            threads.append(t)
            self.set_max_value_and_log(semaphore, limit, f"m{i}", log_values)
            assert semaphore._value <= semaphore.max_value
            assert semaphore.max_value == limit

        for ts in threads:
            for t in ts:
                t.join()

    def start_in_thread(self, func, delay, *args):
        time.sleep(delay)
        thread = threading.Thread(target=func, args=args)
        thread.start()
        return thread

    @staticmethod
    def acquire_and_log(semaphore, n, key, log_values):
        for _ in range(n):
            semaphore.acquire()
        log_values[key] = datetime.utcnow()

    @staticmethod
    def release_and_log(semaphore, n, key, log_values):
        semaphore.release(n)
        log_values[key] = datetime.utcnow()

    @staticmethod
    def set_max_value_and_log(semaphore, n, key, log_values):
        semaphore.max_value = n
        log_values[key] = datetime.utcnow()
