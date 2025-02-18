import threading


class LimitedSemaphore:
    """
    A semaphore with an adjustable limit on the maximum value
    """

    def __init__(self, max_value: int):
        if max_value < 0:
            raise ValueError(f"max_value must be greater or equal to 0, was {max_value}")

        self._max_value = max_value
        self._value = max_value
        self._setter_waiting = False
        self._set_lock = threading.Lock()
        self._condition = threading.Condition(threading.Lock())

    def acquire(self):
        """
        Acquire a semaphore, decrementing internal counter by one.

        If internal counter is greater than zero on entry, return immediately,
        otherwise block and wait until another thread has called release.
        """
        with self._condition:
            # Give priority to the max_value setter and then wait until a semaphore is available
            while self._setter_waiting or self._value == 0:
                self._condition.wait()
            self._value -= 1

    __enter__ = acquire

    def release(self, n: int = 1):
        """
        Release a semaphore, incrementing the internal counter by one or more.

        When the counter goes above the configured maximum value, raises a ValueError.

        When the counter is zero on entry and another thread is waiting for it
        to become larger than zero again, wake up that thread.
        """
        if n < 1:
            raise ValueError(f"n must be one or more, was {n}")

        with self._condition:
            # Increase counter if within limit and notify waiting threads
            if self._value + n > self._max_value:
                raise ValueError(f"LimitedSemaphore released too many times, maximum is set to {self._max_value}")
            self._value += n
            self._condition.notify_all()

    def __exit__(self, t, v, tb):
        self.release()

    @property
    def max_value(self) -> int:
        """
        Maximum value for the internal counter.

        If the semaphore exceeds this a ValueError will be raised.

        Setting the maximum value higher than the current maximum will immediately allow further acquire calls.
        Conversely, reducing the maximum will block until enough threads have released to satisfy the new maximum.
        """
        with self._condition:
            return self._max_value

    @max_value.setter
    def max_value(self, max_value: int):
        with self._set_lock:
            with self._condition:
                diff = max_value - self._max_value
                if diff == 0:
                    return

                # Block new acquires
                self._setter_waiting = True

                # If decreasing max value, wait for enough releases to safely increase max
                # If increasing max value, can just go ahead
                if diff < 0:
                    while diff + self._value < 0:
                        self._condition.wait()

                # Modify the (max) value and notify waiting threads
                self._value += diff
                self._max_value = max_value
                self._setter_waiting = False
                self._condition.notify_all()
