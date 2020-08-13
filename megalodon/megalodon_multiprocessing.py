from time import sleep
import multiprocessing as mp
from collections import namedtuple
from multiprocessing.connection import wait
from multiprocessing.queues import Queue as mpQueue

from megalodon import logging, megalodon_helper as mh


_FULL_SLEEP_TIME = 1

GETTER_PROC = namedtuple('getter_proc', ('queue', 'proc', 'conn'))

LOGGER = logging.get_logger()


###########################
# Multi-processing Helper #
###########################

class CountingMPQueue(mpQueue):
    """ Minimal version of multiprocessing queue maintaining a queue size
    counter
    """

    def __init__(self, **kwargs):
        super().__init__(ctx=mp.get_context(), **kwargs)
        self._size = mp.Value('i', 0)
        self.maxsize = None
        if 'maxsize' in kwargs:
            self.maxsize = kwargs['maxsize']

    def put(self, *args, **kwargs):
        super().put(*args, **kwargs)
        with self._size.get_lock():
            self._size.value += 1

    def get(self, *args, **kwargs):
        rval = super().get(*args, **kwargs)
        with self._size.get_lock():
            self._size.value -= 1
        return rval

    def qsize(self):
        qsize = max(0, self._size.value)
        if self.maxsize is not None:
            return min(self.maxsize, qsize)
        return qsize

    def empty(self):
        return self.qsize() <= 0


def create_getter_q(getter_func, args, max_size=mh._MAX_QUEUE_SIZE):
    if max_size is None:
        q = CountingMPQueue()
    else:
        q = CountingMPQueue(maxsize=max_size)
    main_conn, conn = mp.Pipe()
    p = mp.Process(target=getter_func, daemon=True, args=(q, conn, *args))
    p.start()
    return GETTER_PROC(q, p, main_conn)


class ConnWithSize:
    def __init__(
            self, conn, size, max_size=mh._MAX_QUEUE_SIZE, name='ConnWithSize',
            full_sleep_time=_FULL_SLEEP_TIME):
        if not isinstance(conn, mp.connection.Connection):
            raise mh.MegaError((
                'ConnWithSize initialized with non-connection object. ' +
                'Object type: {}').format(type(conn)))
        if not isinstance(size, mp.sharedctypes.Synchronized) and \
           isinstance(size.value, int):
            raise mh.MegaError((
                'ConnWithSize initialized with non-synchronized size ' +
                'object. Object type: {}').format(type(size)))
        self._conn = conn
        self._size = size
        self.max_size = max_size
        self.full_sleep_time = full_sleep_time
        self.name = name

    def qsize(self):
        return max(min(self._size.value, mh._MAX_QUEUE_SIZE), 0)

    def full(self):
        if self.max_size is None:
            return False
        return self.qsize() >= self.max_size

    def put(self, value):
        # enforce artificial queue max size with dulplex pipes
        if self.full():
            LOGGER.debug('ThrottlingSimplexQueue')
            sleep(self.full_sleep_time)
        with self._size.get_lock():
            self._size.value += 1
        self._conn.send(value)

    def close(self):
        self._conn.close()
        del self._conn


class SimplexManyToOneQueue:
    """ This object is a more efficient version of a multiprocessing.Queue for
    use when many connections will send information in one direction to a
    single connection.

    The get_conn class function will return a ConnWithSize object which can
    send information to be recieved by the wait_recv function of this class.
    """

    def __init__(
            self, return_conns=True, max_size=mh._MAX_QUEUE_SIZE,
            name='SimplexQueue'):
        self.return_conns = return_conns
        self._conns = []
        self._size = mp.Value('i', 0)
        self.max_size = max_size
        self.name = name

    def get_conn(self):
        if not self.return_conns:
            return
        _my_conn, r_conn = mp.Pipe(duplex=False)
        self._conns.append(_my_conn)
        return ConnWithSize(r_conn, self._size, self.max_size, self.name)

    def qsize(self):
        return max(min(self._size.value, mh._MAX_QUEUE_SIZE), 0)

    def empty(self):
        return self.qsize() <= 0

    @property
    def has_valid_conns(self):
        return len(self._conns) > 0

    def wait_recv(self):
        for conn in wait(self._conns):
            try:
                r_val = conn.recv()
                with self._size.get_lock():
                    self._size.value -= 1
            except EOFError:
                # when connection is closed in worker process EOFError is
                # triggered, so remove that connection
                self._conns.remove(conn)
            else:
                yield r_val
